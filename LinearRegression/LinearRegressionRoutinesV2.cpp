#include <iostream>
#include <mpi.h>
#include <limits.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <errno.h>
#include <fstream>

#define DONE_TAG 42424242

#ifndef GANTT
#define GANTT 0
#endif

#ifndef CANCEL_LG
#define CANCEL_LG 0
#endif

#define CANCEL_LG_TAG DONE_TAG+1

using namespace std;

static void get_incomplete_worker(int *worker_iterations, int current_iteration, int n, int *r) {
    //Worker is incomplete if on current_iteration.
    for (int i = 0; i < n; i++) {
	if (worker_iterations[i] != current_iteration+1) {
	    *r = i;
	    break;
	}
    }
}

static void get_both_incomplete_workers(int *worker_iterations, int current_iteration, int n, int *smaller, int *larger) {
    int i;
    for (i = 0; i < n; i++) {
	if (worker_iterations[i] != current_iteration+1) {
	    *smaller = i++;
	    break;
	}
    }
    for (; i < n; i++) {
	if (worker_iterations[i] != current_iteration+1) {
	    *larger = i++;
	    return;
	}
    }
}

//Decode full output and put into output
void decode_parities(int n_completed, double *out, int length, int length_per_worker, double *parity1, double *parity2,
		     double *rest_sum1, double *rest_sum2, int master_iteration, int n_workers, int *worker_iterations) {
    //None missing
    if (n_completed == n_workers) return;

    //1 missing worker
    if (n_completed == n_workers - 1) {
	int worker_incomplete = -1;
	get_incomplete_worker(worker_iterations, master_iteration, n_workers, &worker_incomplete);

	//If the second to last process completed, and all but 1 process completed,
	//then compute the result of the last process
	if (worker_iterations[n_workers-2] > master_iteration) {
	    if (worker_incomplete != n_workers - 1) {
		double * out_i = &out[length_per_worker*worker_incomplete];
		memcpy(out_i, rest_sum1, sizeof(double) * length_per_worker);
	    }
	}
    }
    else if (n_completed == n_workers - 2) {
	int smaller_incomplete = -1, larger_incomplete = -1;
	get_both_incomplete_workers(worker_iterations, master_iteration, n_workers, &smaller_incomplete, &larger_incomplete);
	if (larger_incomplete == n_workers-2) {
	    double * out_i = &out[length_per_worker*smaller_incomplete];
	    vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(smaller_incomplete+1), length_per_worker);
	    memcpy(out_i, rest_sum2, sizeof(double) * length_per_worker);
	}
	if (larger_incomplete == n_workers-1 && smaller_incomplete < n_workers-2) {
	    double * out_i = &out[length_per_worker*smaller_incomplete];
	    memcpy(out_i, rest_sum1, sizeof(double) * length_per_worker);
	}
	if (larger_incomplete < n_workers-2) {
	    double * smaller_out = &out[length_per_worker*smaller_incomplete];
	    double * larger_out = &out[length_per_worker*larger_incomplete];
	    vector_vector_add_scalark(rest_sum2, rest_sum1, rest_sum2, -(smaller_incomplete+1), length_per_worker);
	    vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(larger_incomplete-smaller_incomplete), length_per_worker);
	    vector_vector_subtract(rest_sum1, rest_sum2, rest_sum1, length_per_worker);
	    memcpy(larger_out, rest_sum2, sizeof(double) * length_per_worker);
	    memcpy(smaller_out, rest_sum1, sizeof(double) * length_per_worker);
	}
    }
}

double * master_process_linear_regression(int m_cols, int n_rows, string data_path, double learning_rate) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);
    int m_cols_per_worker = m_cols / (n_workers-2);

    //Master data
    double *out_transpose, *gradient, *out;
    double *labels;
    double *rest_sum1, *rest_sum2, *parity1, *parity2;
    double *dummy_recv;
    int max_row_col = max(n_rows_per_worker, m_cols_per_worker);

    //allocate/load memory
    gradient = (double *)calloc(m_cols, sizeof(double));
    out_transpose = (double *)calloc(n_rows, sizeof(double));
    out = (double *)calloc(m_cols, sizeof(double));
    labels = load_vector(n_rows, "labels.dat");

    rest_sum1 = (double *)calloc(max_row_col, sizeof(double));
    rest_sum2 = (double *)calloc(max_row_col, sizeof(double));
    parity1 = (double *)calloc(max_row_col, sizeof(double));
    parity2 = (double *)calloc(max_row_col, sizeof(double));
    dummy_recv = (double *)malloc(sizeof(double) * max_row_col);

    //Start and end times for workers
    long long int start_times[n_procs][N_ITERS*2+1];
    long long int end_times[n_procs][N_ITERS*2+1];

    //Worker iterations. Note: Iteration refers to a single matrix multiplication
    //not a single iteration of linear regression. This is to simplify
    //the code. So iteration is even if on Ax^(t), and odd if doing A^T(Ax^(t)-y).
    int worker_iterations[n_workers];
    MPI_Request worker_requests[n_workers];
    int master_iteration = 0;
    int n_completed_this_iteration = 0;
    int need_to_send_data = 1;
    for (int i = 0; i < n_workers; i++) worker_iterations[i] = 0;

    //Wait for workers to finish and continue iterations
    while (master_iteration < N_ITERS * 2) {
	//cout << "MASTER ITER: " << master_iteration << " COMPLETED: " << n_completed_this_iteration << endl;
	//cout << "SENDING DATA..." << endl;
	///////////////////////
	// Step 1. Send Data //
	///////////////////////
	if (need_to_send_data) {
	    start_times[0][master_iteration] = get_time();
	    //On Ax^(t)
	    if (master_iteration % 2 == 0) {
		for (int i = 0; i < n_workers; i++) {
		    start_times[i+1][master_iteration] = get_time();
		    MPI_Isend(out, m_cols, MPI_DOUBLE, i+1, master_iteration, MPI_COMM_WORLD, &worker_requests[i]);
		}
	    }
	    //On A^T(Ax^(t)-y)
	    else {
		for (int i = 0; i < n_workers; i++) {
		    start_times[i+1][master_iteration] = get_time();
		    MPI_Isend(out_transpose, n_rows, MPI_DOUBLE, i+1, master_iteration, MPI_COMM_WORLD, &worker_requests[i]);
		}
	    }
	    need_to_send_data = 0;
	}

	//cout << "PROBING AND COLLECTING DATA..." << endl;
	////////////////////////////////////
	// Step 2. Probe and collect data //
	////////////////////////////////////
	for (int i = 0; i < n_workers; i++) {

	    //Check if worker sent data
	    int result_in_queue = 0;
	    MPI_Status message_info;
	    MPI_Iprobe(i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &result_in_queue, &message_info);

	    //Worker sent data
	    if (result_in_queue) {
		if (GANTT && worker_iterations[i] != 0) {
		    end_times[i+1][worker_iterations[i]] = get_time();
		    ofstream gantt_out(("gantt/gantt_" + int_to_string(i+1)).c_str(),  std::fstream::app);
		    gantt_out << start_times[i+1][worker_iterations[i]] << " " << end_times[i+1][worker_iterations[i]] << endl;
		}

		//If worker is in previous iteration, ignore it
		if (worker_iterations[i] < master_iteration) {
		    //recv and move on
		    int count = -1;
		    if (worker_iterations[i] % 2 == 0)
			count = n_rows_per_worker;
		    else
			count = m_cols_per_worker;
		    MPI_Recv(dummy_recv, count, MPI_DOUBLE, i+1, message_info.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    worker_iterations[i]++;
		}
		//If worker is in current iteration, receive values
		else if (worker_iterations[i] == master_iteration) {

		    //These are parity workers
		    if (i+1 >= n_procs-2) {
			int n_vals_per_worker = 0;
			if (master_iteration % 2 == 0)
			    n_vals_per_worker = n_rows_per_worker;
			else
			    n_vals_per_worker = m_cols_per_worker;
			if (i+1 == n_procs-2) {
			    //cout << "RESTSUM1 " << rest_sum1 << endl;
			    MPI_Recv(parity1, n_vals_per_worker, MPI_DOUBLE, i+1, message_info.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			    //cout << "RESTSUM1 " << rest_sum1 << endl;
			    vector_vector_sum(rest_sum1, parity1, rest_sum1, n_vals_per_worker);
			}
			else {
			    MPI_Recv(parity2, n_vals_per_worker, MPI_DOUBLE, i+1, message_info.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			    vector_vector_sum(rest_sum2, parity2, rest_sum2, n_vals_per_worker);
			}
		    }
		    //These are regular workers
		    else {
			//out_transpose = Ax^(t)
			if (master_iteration % 2 == 0) {
			    MPI_Recv(&out_transpose[i*n_rows_per_worker], n_rows_per_worker, MPI_DOUBLE, i+1, message_info.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			    //Auxiliary computation to keep track of rest sums for decoding
			    vector_vector_subtract(rest_sum1, &out_transpose[n_rows_per_worker*i], rest_sum1, n_rows_per_worker);
			    vector_vector_add_scalark(rest_sum2, &out_transpose[n_rows_per_worker*i], rest_sum2, -(i+1), n_rows_per_worker);
			}
			//gradient = A^T(Ax^(t)-y)
			else {
			    MPI_Recv(&gradient[i*m_cols_per_worker], m_cols_per_worker, MPI_DOUBLE, i+1, message_info.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			    //Auxiliary computation to keep track of rest sums for decoding
			    vector_vector_subtract(rest_sum1, &gradient[m_cols_per_worker*i], rest_sum1, m_cols_per_worker);
			    vector_vector_add_scalark(rest_sum2, &gradient[m_cols_per_worker*i], rest_sum2, -(i+1), m_cols_per_worker);
			}
		    }

		    worker_iterations[i]++;
		    n_completed_this_iteration++;
		}
		//If worker is in next iteration, something should be wrong
		else if (worker_iterations[i] > master_iteration) {
		    cout << "Worker is in next iteration of master... Something went wrong? " << "WORKER : " << i << " " << worker_iterations[i] << " vs " << master_iteration << " ACTUAL TAG: " << message_info.MPI_TAG << endl;
		}
	    }
	}

	//cout << "DECODING DATA... " << n_completed_this_iteration << endl;
	///////////////////////////////////////////////////////////////////////////
	// Step 3. Decode partial results using codes. Advance to next iteration //
	///////////////////////////////////////////////////////////////////////////
	if (n_completed_this_iteration >= n_workers-2) {
	    //cout << "ACTUALLY DECODING... " << endl;
	    //Decode
	    if (master_iteration % 2 == 0) {
		decode_parities(n_completed_this_iteration, out_transpose, n_rows, n_rows_per_worker, parity1, parity2,
				rest_sum1, rest_sum2, master_iteration, n_workers, worker_iterations);
	    }
	    else {
		decode_parities(n_completed_this_iteration, gradient, m_cols, m_cols_per_worker, parity1, parity2,
				rest_sum1, rest_sum2, master_iteration, n_workers, worker_iterations);
	    }

	    //cout << "COMPUTATION TO MOVE TO NEXT ITERATION..." << endl;
	    //Do computation to move to next iteration
	    if (master_iteration % 2 == 0) {
		vector_vector_subtract(out_transpose, labels, out_transpose, n_rows);
		memset(rest_sum1, 0, sizeof(double) * max_row_col);
		memset(rest_sum2, 0, sizeof(double) * max_row_col);
		memset(parity1, 0, sizeof(double) * max_row_col);
		memset(parity2, 0, sizeof(double) * max_row_col);
	    }
	    else {
		//cout << "COMPUTATION MOD 1 " << endl;
		vector_vector_subtract_scalark(out, gradient, out, learning_rate, m_cols);
		//cout << "CLEARING..." << endl;
		memset(out_transpose, 0, sizeof(double) * n_rows);
		memset(gradient, 0, sizeof(double) * m_cols);
		memset(rest_sum1, 0, sizeof(double) * max_row_col);
		memset(rest_sum2, 0, sizeof(double) * max_row_col);
		memset(parity1, 0, sizeof(double) * max_row_col);
		memset(parity2, 0, sizeof(double) * max_row_col);
	    }
	    //cout << "MOVED TO NEXT ITERATION... " << endl;

	    if (GANTT && master_iteration != 0) {
		end_times[0][master_iteration] = get_time();
		ofstream gantt_out(("gantt/gantt_" + int_to_string(0)).c_str(),  std::fstream::app);
		gantt_out << start_times[0][master_iteration] << " " << end_times[0][master_iteration] << endl;
	    }

	    //On next iteration now
	    master_iteration++;
	    n_completed_this_iteration = 0;
	    need_to_send_data = 1;

	    //Cancel those previous iterations
	    if (CANCEL_LG) {
		for (int i = 0; i < n_workers; i++) {
		    //Send an int representing iteration to cancel
		    if (worker_iterations[i] < master_iteration) {
			int cancel_message = worker_iterations[i];
			MPI_Request cancel_request;
			MPI_Isend(&cancel_message, 1, MPI_INT, i+1, CANCEL_LG_TAG, MPI_COMM_WORLD, &cancel_request);
		    }
		}
	    }
	}
    }

    //cout << "DONE! " << endl;

    //Send done message
    for (int i = 0; i < n_workers; i++) {
	int done = 1;
	MPI_Request r;
	MPI_Isend(&done, 1, MPI_INT, i+1, DONE_TAG, MPI_COMM_WORLD, &r);
    }

    //Free mem

    free(gradient);
    free(out_transpose);
    free(labels);
    free(rest_sum1);
    free(rest_sum2);
    free(parity1);
    free(parity2);

    return out;
}

double * worker_process_linear_regression_no_cancel(int m_cols, int n_rows, string data_path, double learning_rate) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);
    int m_cols_per_worker = m_cols / (n_workers-2);

    //Worker process data
    double *sub_data;
    double *sub_vec, *sub_out;
    double *sub_data_transpose;
    double *sub_vec_transpose, *sub_out_transpose;

    //Allocate memory for work
    sub_vec = (double *)calloc(m_cols, sizeof(double));
    sub_out = (double *)calloc(n_rows_per_worker, sizeof(double));
    sub_vec_transpose = (double *)calloc(n_rows, sizeof(double));
    sub_out_transpose = (double *)calloc(m_cols_per_worker, sizeof(double));
    sub_data = load_submatrix(proc_id, m_cols, n_rows_per_worker);
    sub_data_transpose = load_submatrix_t(proc_id, n_rows, m_cols_per_worker);

    //Current worker iteration
    int cur_worker_iteration = 0;

    //Listen to done with whole process
    int done = 0;
    MPI_Request req;
    MPI_Irecv(&done, 1, MPI_INT, 0, DONE_TAG, MPI_COMM_WORLD, &req);

    //While not done yet, listen to messages from master and do computation
    while (!done) {
	int message_exists = 0;
	MPI_Status status;
	MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &message_exists, &status);

	if (message_exists && status.MPI_TAG == cur_worker_iteration) {
	    if (cur_worker_iteration % 2 == 0) {
		MPI_Recv(sub_vec, m_cols, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		matrix_vector_multiply(sub_data, sub_vec, sub_out, m_cols, n_rows_per_worker);
		MPI_Send(sub_out, n_rows_per_worker, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD);
		cur_worker_iteration++;
	    }
	    else {
		MPI_Recv(sub_vec_transpose, n_rows, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		matrix_vector_multiply(sub_data_transpose, sub_vec_transpose, sub_out_transpose, n_rows, m_cols_per_worker);
		MPI_Send(sub_out_transpose, m_cols_per_worker, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD);
		cur_worker_iteration++;

		memset(sub_out, 0, sizeof(double) * n_rows_per_worker);
		memset(sub_out_transpose, 0, sizeof(double) * m_cols_per_worker);
	    }
	}
    }

    //Free memory
    free(sub_data);
    free(sub_vec);
    free(sub_out);
    free(sub_data_transpose);
    free(sub_vec_transpose);
    free(sub_out_transpose);

    return NULL;
}


double * worker_process_linear_regression_cancel(int m_cols, int n_rows, string data_path, double learning_rate) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);
    int m_cols_per_worker = m_cols / (n_workers-2);

    //Worker process data
    double *sub_data;
    double *sub_vec, *sub_out;
    double *sub_data_transpose;
    double *sub_vec_transpose, *sub_out_transpose;

    //Allocate memory for work
    sub_vec = (double *)calloc(m_cols, sizeof(double));
    sub_out = (double *)calloc(n_rows_per_worker, sizeof(double));
    sub_vec_transpose = (double *)calloc(n_rows, sizeof(double));
    sub_out_transpose = (double *)calloc(m_cols_per_worker, sizeof(double));
    sub_data = load_submatrix(proc_id, m_cols, n_rows_per_worker);
    sub_data_transpose = load_submatrix_t(proc_id, n_rows, m_cols_per_worker);

    //Current worker iteration
    int cur_worker_iteration = 0;

    //Listen to done with whole process
    int done = 0;
    MPI_Request req;
    MPI_Irecv(&done, 1, MPI_INT, 0, DONE_TAG, MPI_COMM_WORLD, &req);

    //While not done yet, listen to messages from master and do computation
    while (!done) {
	int message_exists = 0;
	MPI_Status status;
	MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &message_exists, &status);

	if (message_exists && status.MPI_TAG == cur_worker_iteration) {

	    //Listen to cancellation signals
	    MPI_Request cancel_req;
	    int cancelled = INT_MAX;
	    MPI_Irecv(&cancelled, 1, MPI_INT, 0, CANCEL_LG_TAG, MPI_COMM_WORLD, &cancel_req);

	    double *work_vec, *data_mat, *out_vec;
	    int n_vals = 0, result_length = 0;

	    if (cur_worker_iteration % 2 == 0) {
		MPI_Recv(sub_vec, m_cols, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		data_mat = sub_data;
		work_vec = sub_vec;
		out_vec = sub_out;
		n_vals = m_cols;
		result_length = n_rows_per_worker;
	    }
	    else {
		MPI_Recv(sub_vec_transpose, n_rows, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		data_mat = sub_data_transpose;
		work_vec = sub_vec_transpose;
		out_vec = sub_out_transpose;
		n_vals = n_rows;
		result_length = m_cols_per_worker;
	    }

	    int computation_done = 0;
	    struct matrix_data_construct data = {proc_id, cur_worker_iteration, data_mat, work_vec, out_vec, n_vals, result_length, &computation_done};
	    pthread_t computation_thread;
	    int pthread_status = pthread_create(&computation_thread, NULL, matrix_computation, (void *)&data);
	    int i = 0, flag = 0;
	    while (!computation_done) {
		usleep(10000);
		if (++i % 2 == 0)
		    MPI_Iprobe(0, CANCEL_LG_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
		//Found a cancel for the current iteration, kill thread and exit
		if (cancelled == cur_worker_iteration) {
		    //cout << "CANCEL_LGLED PROCESS: " << proc_id << " ON ITERATION: " << cancelled << endl;
		    pthread_cancel(computation_thread);
		    break;
	      }
	    }
	    pthread_join(computation_thread, NULL);

	    //matrix_vector_multiply(data_mat, work_vec, out_vec, n_vals, result_length);
	    MPI_Send(out_vec, result_length, MPI_DOUBLE, 0, cur_worker_iteration, MPI_COMM_WORLD);

	    if (cur_worker_iteration % 2 == 1) {
		memset(sub_out, 0, sizeof(double) * n_rows_per_worker);
		memset(sub_out_transpose, 0, sizeof(double) * m_cols_per_worker);
	    }
	    cur_worker_iteration++;

	}
    }

    //Free memory
    free(sub_data);
    free(sub_vec);
    free(sub_out);
    free(sub_data_transpose);
    free(sub_vec_transpose);
    free(sub_out_transpose);

    return NULL;
}


double * worker_process_linear_regression(int m_cols, int n_rows, string data_path, double learning_rate) {
    if (!CANCEL_LG)
	worker_process_linear_regression_no_cancel(m_cols, n_rows, data_path, learning_rate);
    else
	worker_process_linear_regression_cancel(m_cols, n_rows, data_path, learning_rate);
    return NULL;
}


double * coded2_linear_regression_v2(int m_cols, int n_rows, string data_path, double learning_rate) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);
    int m_cols_per_worker = m_cols / (n_workers-2);
    MPI_Barrier(MPI_COMM_WORLD);
    double *ret;
    if (proc_id == 0)
	ret = master_process_linear_regression(m_cols, n_rows, data_path, learning_rate);
    else
	ret = worker_process_linear_regression(m_cols, n_rows, data_path, learning_rate);

    MPI_Barrier(MPI_COMM_WORLD);
    return ret;
}
