#include <iostream>
#include <mpi.h>
#include "../util.cpp"
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <errno.h>
#include <fstream>

#define CANCEL_TAG 320000

#ifndef GANTT
#define GANTT 0
#endif

struct matrix_data_construct {
    int proc_id, iteration;
    double *submat, *subvec, *subout;
    int m_cols, n_rows;
    int *completed;
};

void get_incomplete_worker(int *completed, int n_workers, int *out) {
    for (int i = 0; i < n_workers; i++) {
	if (!completed[i]) {
	    *out = i;
	    return;
	}
    }
}

void get_both_incomplete_workers(int *completed, int n_workers, int *smaller, int *larger) {
    int i;
    for (i = 0; i < n_workers; i++) {
	if (!completed[i]) {
	    *smaller = i++;
	    break;
	}
    }
    for (; i < n_workers; i++) {
	if (!completed[i]) {
	    *larger = i;
	    return;
	}
    }
}

void extreme_matrix_vector_multiply_parallel(int m_cols, int n_rows) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows;

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout, *parity;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Main work
    if (proc_id == 0) {
	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows);
	long long int times[n_workers];

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    times[i-1] = get_time();
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	}
	MPI_Request_free(&send_req);

	//Collect results
	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];

	//Collect sub results
	for (int i = 1; i < n_procs; i++) {
	    MPI_Irecv(out, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	}

	int completed_worker = -1;
	MPI_Waitany(n_workers, recv_reqs, &completed_worker, MPI_STATUS_IGNORE);

	long long int t2 = get_time();
	check_correct(out, n_rows);
	cout << "-------------------------------------------" << endl;
	cout << "EXTREME MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	int n_done = 0, completed[n_workers];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_done += 1;
		    }
		}
	    }
	}
	print_worker_roundtrip_times(times, n_workers);
	cout << "-------------------------------------------" << endl;

	free(vec);
	free(out);
    }
    else {
	//Receive the  vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}

void coded2_matrix_vector_multiply_parallel(int m_cols, int n_rows) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout, *parity1, *parity2;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Main work
    if (proc_id == 0) {
	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows);
	double *rest_sum1 = (double *)calloc(n_rows_per_worker, sizeof(double));
	double *rest_sum2 = (double *)calloc(n_rows_per_worker, sizeof(double));
	double *parity1 = (double *)calloc(n_rows_per_worker, sizeof(double));
	double *parity2 = (double *)calloc(n_rows_per_worker, sizeof(double));
	long long int times[n_workers];

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    times[i-1] = get_time();
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	}
	MPI_Request_free(&send_req);

	//Collect results
	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];

	//Collect sub results
	for (int i = 1; i < n_procs; i++) {
	    double * out_i = &out[n_rows_per_worker*(i-1)];
	    if (i == n_procs-2) {
		MPI_Irecv(parity1, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else if (i == n_procs-1) {
		MPI_Irecv(parity2, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else {
		MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	}

	//Wait for n-1 workers to complete, computing partial sums as they arrive
	int n_done = 0, completed[n_workers];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers-2) { //CHANGE BACK
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_done += 1;
			if (i == n_workers-2) {
			    vector_vector_sum(rest_sum1, parity1, rest_sum1, n_rows_per_worker);
			}
			else if (i == n_workers-1) {
			    vector_vector_sum(rest_sum2, parity2, rest_sum2, n_rows_per_worker);
			}
			else{
			    vector_vector_subtract(rest_sum1, &out[n_rows_per_worker*i], rest_sum1, n_rows_per_worker);
			    vector_vector_add_scalark(rest_sum2, &out[n_rows_per_worker*i], rest_sum2, -(i+1), n_rows_per_worker);
			}
		    }
		}
	    }
	}

	//If only 1 missing worker
	if (n_done == n_workers-1) {
	    int worker_incomplete = -1;
	    get_incomplete_worker(completed, n_workers, &worker_incomplete);
	    //If the second to last process completed, and all but 1 process completed,
	    //then compute the result of the last process
	    if (completed[n_workers-2]) {
		if (worker_incomplete != n_workers - 1) {
		    double * out_i = &out[n_rows_per_worker*worker_incomplete];
		    memcpy(out_i, rest_sum1, sizeof(double) * n_rows_per_worker);
		}
	    }
	}
	else if (n_done == n_workers-2) {
	    int smaller_incomplete, larger_incomplete;
	    get_both_incomplete_workers(completed, n_workers, &smaller_incomplete, &larger_incomplete);
	    if (larger_incomplete == n_workers-2) {
		double * out_i = &out[n_rows_per_worker*smaller_incomplete];
		vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(smaller_incomplete+1), n_rows_per_worker);
		memcpy(out_i, rest_sum2, sizeof(double) * n_rows_per_worker);
	    }
	    if (larger_incomplete == n_workers-1 && smaller_incomplete < n_workers-2) {
		double * out_i = &out[n_rows_per_worker*smaller_incomplete];
		memcpy(out_i, rest_sum1, sizeof(double) * n_rows_per_worker);
	    }
	    if (larger_incomplete < n_workers-2) {
		double * smaller_out = &out[n_rows_per_worker*smaller_incomplete];
		double * larger_out = &out[n_rows_per_worker*larger_incomplete];
		vector_vector_add_scalark(rest_sum2, rest_sum1, rest_sum2, -(smaller_incomplete+1), n_rows_per_worker);
		vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(larger_incomplete-smaller_incomplete), n_rows_per_worker);
		vector_vector_subtract(rest_sum1, rest_sum2, rest_sum1, n_rows_per_worker);
		memcpy(larger_out, rest_sum2, sizeof(double) * n_rows_per_worker);
		memcpy(smaller_out, rest_sum1, sizeof(double) * n_rows_per_worker);
	    }
	}

	long long int t2 = get_time();
	check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << "TWO CODED MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_done += 1;
		    }
		}
	    }
	}
	print_worker_roundtrip_times(times, n_workers);
	cout << "-------------------------------------------" << endl;

	free(vec);
	free(out);
	free(rest_sum1);
	free(rest_sum2);
	free(parity1);
	free(parity2);
    }
    else {
	//Receive the  vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}

void coded_matrix_vector_multiply_parallel(int m_cols, int n_rows) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-1);

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout, *parity;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    //Main work
    if (proc_id == 0) {
	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows);
	double *rest_sum = (double *)calloc(n_rows_per_worker, sizeof(double));
	double *parity = (double *)calloc(n_rows_per_worker, sizeof(double));
	long long int times[n_workers];

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    times[i-1] = get_time();
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	}
	MPI_Request_free(&send_req);

	//Collect results
	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];

	//Collect sub results
	for (int i = 1; i < n_procs; i++) {
	    double * out_i = &out[n_rows_per_worker*(i-1)];
	    if (i == n_procs-1) {
		MPI_Irecv(parity, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else {
		MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	}

	//Wait for n-1 workers to complete, computing partial sums as they arrive
	int n_done = 0, completed[n_workers];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers-1) { //CHANGNE BACK
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_done += 1;
			if (i == n_workers-1) {
			    vector_vector_sum(parity, rest_sum, rest_sum, n_rows_per_worker);
			}
			else{
			    vector_vector_subtract(rest_sum, &out[n_rows_per_worker*i], rest_sum, n_rows_per_worker);
			}
		    }
		}
	    }
	}

	//If the last process completed, and not all processes completed,
	//then compute the result of the last process
	if (completed[n_workers-1] && n_done < n_workers) {
	    int worker_incomplete = -1;
	    get_incomplete_worker(completed, n_workers, &worker_incomplete);
	    double * out_i = &out[n_rows_per_worker*worker_incomplete];
	    memcpy(out_i, rest_sum, sizeof(double) * n_rows_per_worker);
	}

	long long int t2 = get_time();
	check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << "ONE CODED MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_done += 1;
		    }
		}
	    }
	}
	print_worker_roundtrip_times(times, n_workers);
	cout << "-------------------------------------------" << endl;

	free(vec);
	free(out);
	free(rest_sum);
	free(parity);
    }
    else {
	//Receive the  vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}

void * matrix_computation(void *data) {
    int type, state;
    //stick_this_thread_to_core(1);
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &type);
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &state);
    struct matrix_data_construct d = *(struct matrix_data_construct *)data;
    *(d.completed) = 0;
    if ((d.proc_id == 3) && d.iteration == 0) {
	sleep(5);
    }
    matrix_vector_multiply(d.submat, d.subvec, d.subout, d.m_cols, d.n_rows);
    *(d.completed) = 1;
    return NULL;
}

//Same as coded2_matrix_vector_multiply_parallel, but does not wait for the last
//2 processes to stop after getting enough results to compute the answer
void coded2_matrix_vector_multiply_parallel_singlerun(int m_cols, int n_rows, int tag, double *submat, double *subvec, double *subout,
						      double *vec, double *out, double *outkk, double *rest_sum1, double *rest_sum2, double *parity1, double *parity2) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);
    int cancelled = -1;
    long long int start_wrk_time, end_wrk_time;
    MPI_Request recv;

    if (proc_id != 0) {
      //Check to see if this process was cancelled
      MPI_Irecv(&cancelled, 1, MPI_INT, 0, CANCEL_TAG, MPI_COMM_WORLD, &recv);
    }

    if (tag == 0)
      MPI_Barrier(MPI_COMM_WORLD);
    if (GANTT)
	    start_wrk_time = get_time();

    //Main work
    if (proc_id == 0) {
	//load_input_vector_no_alloc(m_cols, vec);

	//Distribute vectors
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &send_req);
	}
	MPI_Request_free(&send_req);

	//Collect results
	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];

	//Collect sub results
	for (int i = 1; i < n_procs; i++) {
	    double * out_i = &outkk[n_rows_per_worker*(i-1)];
	    if (i == n_procs-2) {
		MPI_Irecv(parity1, n_rows_per_worker, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else if (i == n_procs-1) {
		MPI_Irecv(parity2, n_rows_per_worker, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else {
		MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	}

	//Wait for n-1 workers to complete, computing partial sums as they arrive
	int n_done = 0, completed[n_workers];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers-2) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
		      completed[i] = 1;
			n_done += 1;
			if (i == n_workers-2) {
			    vector_vector_sum(rest_sum1, parity1, rest_sum1, n_rows_per_worker);
			}
			else if (i == n_workers-1) {
			    vector_vector_sum(rest_sum2, parity2, rest_sum2, n_rows_per_worker);
			}
			else{
			    memcpy(&out[n_rows_per_worker*i], &outkk[n_rows_per_worker*i], sizeof(double)*n_rows_per_worker);
			    vector_vector_subtract(rest_sum1, &out[n_rows_per_worker*i], rest_sum1, n_rows_per_worker);
			    vector_vector_add_scalark(rest_sum2, &out[n_rows_per_worker*i], rest_sum2, -(i+1), n_rows_per_worker);
			}
		    }
		}
	    }
	}

	//If only 1 missing worker
	if (n_done == n_workers-1) {
	    int worker_incomplete = -1;
	    get_incomplete_worker(completed, n_workers, &worker_incomplete);

	    //If the second to last process completed, and all but 1 process completed,
	    //then compute the result of the last process
	    if (completed[n_workers-2]) {
		if (worker_incomplete != n_workers - 1) {
		    double * out_i = &out[n_rows_per_worker*worker_incomplete];
		    memcpy(out_i, rest_sum1, sizeof(double) * n_rows_per_worker);
		}
	    }
	}
	else if (n_done == n_workers-2) {
	    int smaller_incomplete, larger_incomplete;
	    get_both_incomplete_workers(completed, n_workers, &smaller_incomplete, &larger_incomplete);
	    if (larger_incomplete == n_workers-2) {
		double * out_i = &out[n_rows_per_worker*smaller_incomplete];
		vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(smaller_incomplete+1), n_rows_per_worker);
		memcpy(out_i, rest_sum2, sizeof(double) * n_rows_per_worker);
	    }
	    if (larger_incomplete == n_workers-1 && smaller_incomplete < n_workers-2) {
		double * out_i = &out[n_rows_per_worker*smaller_incomplete];
		memcpy(out_i, rest_sum1, sizeof(double) * n_rows_per_worker);
	    }
	    if (larger_incomplete < n_workers-2) {
		double * smaller_out = &out[n_rows_per_worker*smaller_incomplete];
		double * larger_out = &out[n_rows_per_worker*larger_incomplete];
		vector_vector_add_scalark(rest_sum2, rest_sum1, rest_sum2, -(smaller_incomplete+1), n_rows_per_worker);
		vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(larger_incomplete-smaller_incomplete), n_rows_per_worker);
		vector_vector_subtract(rest_sum1, rest_sum2, rest_sum1, n_rows_per_worker);
		memcpy(larger_out, rest_sum2, sizeof(double) * n_rows_per_worker);
		memcpy(smaller_out, rest_sum1, sizeof(double) * n_rows_per_worker);
	    }
	}

	//Go through and cancel requests.
	for (int i = 0; i < n_workers; i++) {
	  if (!completed[i]) {
	      int did_complete_now = 0;
	      MPI_Test(&recv_reqs[i], &did_complete_now, MPI_STATUS_IGNORE);
	      if (!did_complete_now) {
		  //Send to that worker that we are cancelling his request
		MPI_Request cancel_request;

		//MPI_Rsend(&cancel, 1, MPI_INT, i+1, CANCEL_TAG, MPI_COMM_WORLD);
		//cout << "CANCEL " << i+1 << " iter " << tag << endl;
		MPI_Isend(&tag, 1, MPI_INT, i+1, CANCEL_TAG, MPI_COMM_WORLD, &cancel_request);

		//Cancel
		//MPI_Cancel(&recv_reqs[i]);
	      }
	  }
	}

	if (GANTT) {
	    end_wrk_time = get_time();
	    ofstream gantt_out(("gantt/gantt_" + int_to_string(proc_id)).c_str(),  std::fstream::app);
	    gantt_out << start_wrk_time << " " << end_wrk_time << endl;
	}
    }
    else {

	//Receive the vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	//Pack the data into a single construct
	int computation_completed = 0;
	struct matrix_data_construct data = {proc_id, tag, submat, subvec, subout, m_cols, n_rows_per_worker, &computation_completed};

	//Perform computation on separate thread
	pthread_t computation_thread;
	int pthread_status = pthread_create(&computation_thread, NULL, matrix_computation, (void *)&data);
	//Wait for computation thread to finish, or for it to be cancelled
	int i = 0;
	int flag = -1;
	//pthread_detach(computation_thread);
	//cout << "COMPUTATION " << computation_completed << endl;
	/*while (computation_completed == 0) {
	  usleep(10000);
	  if (++i % 2 == 0)
	    //MPI_Test(&recv, &flag, MPI_STATUS_IGNORE);
	    MPI_Iprobe(0, CANCEL_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
	  //Found a cancel for the current iteration, kill thread and exit
	  if (cancelled == tag) {
	    //cout << "CANCELLED PROCESS: " << proc_id << " ON ITERATION: " << cancelled << endl;
	    pthread_cancel(computation_thread);
	    break;
	  }
	  }*/
	pthread_join(computation_thread, NULL);

	if (GANTT) {
	    end_wrk_time = get_time();
	    ofstream gantt_out(("gantt/gantt_" + int_to_string(proc_id)).c_str(),  std::fstream::app);
	    gantt_out << start_wrk_time << " " << end_wrk_time << endl;
	}

	//Send the result
	if (cancelled != tag) {
	    MPI_Request req;
	    MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	    MPI_Cancel(&recv);
	}

	//Wait for response
	//MPI_Wait(&recv, MPI_STATUS_IGNORE);
	//if (cancelled != tag);
    }
}

void coded2_matrix_vector_multiply_parallel_multirun(int m_cols, int n_rows, int num_runs) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);

    double *vec, *out, *outkk, *rest_sum1, *rest_sum2, *parity1, *parity2;
    double *submat, *subvec, *subout;

    if (proc_id == 0) {
	out = generate_empty_vector(n_rows);
	outkk = generate_empty_vector(n_rows);
	rest_sum1 = (double *)calloc(n_rows_per_worker, sizeof(double));
	rest_sum2 = (double *)calloc(n_rows_per_worker, sizeof(double));
	parity1 = (double *)calloc(n_rows_per_worker, sizeof(double));
	parity2 = (double *)calloc(n_rows_per_worker, sizeof(double));
	vec = generate_empty_vector(m_cols);
    }
    else {
	submat = (double *)malloc(sizeof(double) * m_cols * n_rows);
	load_submatrix_no_alloc(proc_id, m_cols, n_rows_per_worker, submat);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    long long int t1 = get_time();
    for (int i = 0; i < num_runs; i++) {
	coded2_matrix_vector_multiply_parallel_singlerun(m_cols, n_rows, i,
							 submat, subvec, subout,
							 vec, out, outkk, rest_sum1, rest_sum2,
							 parity1, parity2);
	if (proc_id == 0) {
	    memset(out, 0, sizeof(double) * n_rows);
	    memset(rest_sum1, 0, sizeof(double) * n_rows_per_worker);
	    memset(rest_sum2, 0, sizeof(double) * n_rows_per_worker);
	    memset(parity1, 0, sizeof(double) * n_rows_per_worker);
	    memset(parity2, 0, sizeof(double) * n_rows_per_worker);
	}
	else {
	    memset(subout, 0, sizeof(double) * n_rows_per_worker);
	}
    }
    long long int t2 = get_time();
    if (proc_id == 0)
	cout << "ELAPSED TIME: " << t2-t1 << endl;
}
