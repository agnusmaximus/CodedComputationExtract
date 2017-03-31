//Same methods as in coded_matrix_vector_multiply, but profiles
//the processing and communication times of workers, etc

#include <iostream>
#include <mpi.h>
#include "../util.cpp"
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <string.h>

void coded_matrix_vector_multiply_parallel_profile(int m_cols, int n_rows) {

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

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    //The last element of the input vector is the time
	    vec[m_cols-1] = (double)get_time();
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
	double worker_stats[n_workers][3];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers-1) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			completed[i] = 1;
			n_done += 1;
			if (i == n_workers-1) {
			    worker_stats[i][0] = (double)get_time() - parity[n_rows_per_worker-1]; //Time to send solution
			    worker_stats[i][1] = parity[n_rows_per_worker-2];
			    worker_stats[i][2] = parity[n_rows_per_worker-3];
			    vector_vector_sum(parity, rest_sum, rest_sum, n_rows_per_worker);
			}
			else{
			    worker_stats[i][0] = (double)get_time() - out[n_rows_per_worker*(i+1)-1];
			    worker_stats[i][1] = out[n_rows_per_worker*(i+1)-2];
			    worker_stats[i][2] = out[n_rows_per_worker*(i+1)-3];
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
	    //Do the work, but do not write to output; wait for worker to finish
	    //memcpy(out_i, rest_sum, sizeof(double) * n_rows_per_worker);
	}

	//Don't check correctness...
	//check_correct(out, n_rows);

	long long int t2 = get_time();
	cout << "-------------------------------------------" << endl;
	cout << "ONE CODED MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	int late_worker = -1;
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			completed[i] = 1;
			n_done += 1;
			late_worker = i+1;
			if (i == n_workers-1) {
			    worker_stats[i][0] = (double)get_time() - parity[n_rows_per_worker-1]; //Time to send solution
			    worker_stats[i][1] = parity[n_rows_per_worker-2];
			    worker_stats[i][2] = parity[n_rows_per_worker-3];
			}
			else {
			    worker_stats[i][0] = (double)get_time() - out[n_rows_per_worker*(i+1)-1]; //Time to send solution
			    worker_stats[i][1] = out[n_rows_per_worker*(i+1)-2];
			    worker_stats[i][2] = out[n_rows_per_worker*(i+1)-3];
			}
		    }
		}
	    }
	}

	//Print out worker stats
	for (int i = 1; i < n_procs; i++) {
	    double time_send_solution = worker_stats[i-1][0];
	    double time_receive_input = worker_stats[i-1][1];
	    double work_time = worker_stats[i-1][2];
	    if (i == late_worker) cout << "LATE ";
	    cout << "WORKER " << i << endl << " RECEIVE INPUT TIME: " << time_receive_input << endl << " WORK TIME: " << work_time << endl << " SEND RESULT TIME: " << time_send_solution << endl;
	}
	cout << "-------------------------------------------" << endl;

	free(vec);
	free(out);
	free(rest_sum);
	free(parity);
    }
    else {
	//Receive the  vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//The last element of the input vector is the time that it was sent
	long long int time_sent = (long long int)subvec[m_cols-1];
	long long int time_received = get_time();
	long long int time_to_send_to_worker = time_received - time_sent;
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	long long int time_to_compute = get_time() - time_received;
	long long int time_resent = get_time();
	//Last element in returned vector is the time to send the solution
	subout[n_rows_per_worker-1] = (double)time_resent;
	//Second to last element in returned vector is the time to receive the input
	subout[n_rows_per_worker-2] = (double)time_to_send_to_worker;
	//Third to last element in returned vector is the time for the worker to do the computation
	subout[n_rows_per_worker-3] = (double)time_to_compute;
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}

void coded2_matrix_vector_multiply_parallel_profile(int m_cols, int n_rows) {

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

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    vec[m_cols-1] = (double)get_time();
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
	double worker_stats[n_workers][3];
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
			    worker_stats[i][0] = (double)get_time() - parity1[n_rows_per_worker-1]; //Time to send solution
			    worker_stats[i][1] = parity1[n_rows_per_worker-2];
			    worker_stats[i][2] = parity1[n_rows_per_worker-3];
			    vector_vector_sum(rest_sum1, parity1, rest_sum1, n_rows_per_worker);
			}
			else if (i == n_workers-1) {
			    worker_stats[i][0] = (double)get_time() - parity2[n_rows_per_worker-1]; //Time to send solution
			    worker_stats[i][1] = parity2[n_rows_per_worker-2];
			    worker_stats[i][2] = parity2[n_rows_per_worker-3];
			    vector_vector_sum(rest_sum2, parity2, rest_sum2, n_rows_per_worker);
			}
			else{
			    worker_stats[i][0] = (double)get_time() - out[n_rows_per_worker*(i+1)-1];
			    worker_stats[i][1] = out[n_rows_per_worker*(i+1)-2];
			    worker_stats[i][2] = out[n_rows_per_worker*(i+1)-3];
			    vector_vector_subtract(rest_sum1, &out[n_rows_per_worker*i], rest_sum1, n_rows_per_worker);
			    vector_vector_add_scalark(rest_sum2, &out[n_rows_per_worker*i], rest_sum2, -(i+1), n_rows_per_worker);
			}
		    }
		}
	    }
	}

	int late_worker1 = -1, late_worker2 = -1;

	//If only 1 missing worker
	if (n_done == n_workers-1) {
	    int worker_incomplete = -1;
	    get_incomplete_worker(completed, n_workers, &worker_incomplete);
	    late_worker1 = worker_incomplete+1;
	    //If the second to last process completed, and all but 1 process completed,
	    //then compute the result of the last process
	    if (completed[n_workers-2]) {
		if (worker_incomplete != n_workers - 1) {
		    double * out_i = &out[n_rows_per_worker*worker_incomplete];

		    //memcpy(out_i, rest_sum1, sizeof(double) * n_rows_per_worker);
		}
	    }
	}
	else if (n_done == n_workers-2) {
	    int smaller_incomplete, larger_incomplete;
	    get_both_incomplete_workers(completed, n_workers, &smaller_incomplete, &larger_incomplete);
	    late_worker1 = smaller_incomplete+1;
	    late_worker2 = larger_incomplete+1;
	    if (larger_incomplete == n_workers-2) {
		double * out_i = &out[n_rows_per_worker*smaller_incomplete];
		vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(smaller_incomplete+1), n_rows_per_worker);
		//memcpy(out_i, rest_sum2, sizeof(double) * n_rows_per_worker);
	    }
	    if (larger_incomplete == n_workers-1 && smaller_incomplete < n_workers-2) {
		double * out_i = &out[n_rows_per_worker*smaller_incomplete];
		//memcpy(out_i, rest_sum1, sizeof(double) * n_rows_per_worker);
	    }
	    if (larger_incomplete < n_workers-2) {
		double * smaller_out = &out[n_rows_per_worker*smaller_incomplete];
		double * larger_out = &out[n_rows_per_worker*larger_incomplete];
		vector_vector_add_scalark(rest_sum2, rest_sum1, rest_sum2, -(smaller_incomplete+1), n_rows_per_worker);
		vector_vector_scalar_multiply(rest_sum2, rest_sum2, 1/(double)(larger_incomplete-smaller_incomplete), n_rows_per_worker);
		vector_vector_subtract(rest_sum1, rest_sum2, rest_sum1, n_rows_per_worker);
		//memcpy(larger_out, rest_sum2, sizeof(double) * n_rows_per_worker);
		//memcpy(smaller_out, rest_sum1, sizeof(double) * n_rows_per_worker);
	    }
	}

	//check_correct(out, n_rows);

	long long int t2 = get_time();
	cout << "-------------------------------------------" << endl;
	cout << "TWO CODED MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			completed[i] = 1;
			n_done += 1;
			if (i == n_workers-2) {
			    worker_stats[i][0] = (double)get_time() - parity1[n_rows_per_worker-1]; //Time to send solution
			    worker_stats[i][1] = parity1[n_rows_per_worker-2];
			    worker_stats[i][2] = parity1[n_rows_per_worker-3];
			}
			else if (i == n_workers-1) {
			    worker_stats[i][0] = (double)get_time() - parity2[n_rows_per_worker-1]; //Time to send solution
			    worker_stats[i][1] = parity2[n_rows_per_worker-2];
			    worker_stats[i][2] = parity2[n_rows_per_worker-3];
			}
			else{
			    worker_stats[i][0] = (double)get_time() - out[n_rows_per_worker*(i+1)-1];
			    worker_stats[i][1] = out[n_rows_per_worker*(i+1)-2];
			    worker_stats[i][2] = out[n_rows_per_worker*(i+1)-3];
			}
		    }
		}
	    }
	}
	//Print out worker stats
	for (int i = 1; i < n_procs; i++) {
	    double time_send_solution = worker_stats[i-1][0];
	    double time_receive_input = worker_stats[i-1][1];
	    double work_time = worker_stats[i-1][2];
	    if (i == late_worker1 || i == late_worker2) cout << "LATE ";
	    cout << "WORKER " << i << endl << " RECEIVE INPUT TIME: " << time_receive_input << endl << " WORK TIME: " << work_time << endl << " SEND RESULT TIME: " << time_send_solution << endl;
	}
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
	//The last element of the input vector is the time that it was sent
	long long int time_sent = (long long int)subvec[m_cols-1];
	long long int time_received = get_time();
	long long int time_to_send_to_worker = time_received - time_sent;
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	long long int time_to_compute = get_time() - time_received;
	long long int time_resent = get_time();
	//Last element in returned vector is the time to send the solution
	subout[n_rows_per_worker-1] = (double)time_resent;
	//Second to last element in returned vector is the time to receive the input
	subout[n_rows_per_worker-2] = (double)time_to_send_to_worker;
	//Third to last element in returned vector is the time for the worker to do the computation
	subout[n_rows_per_worker-3] = (double)time_to_compute;
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}

void uncoded_matrix_vector_multiply_parallel_profile(int m_cols, int n_rows) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / n_workers;

    if (proc_id == n_procs-1)
	n_rows_per_worker = n_rows-n_rows_per_worker*(proc_id-1);

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Main parallel code
    if (proc_id == 0) {

	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows);

	//Distribute vectors
	int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    vec[m_cols-1] = (double)get_time();
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	}
	MPI_Request_free(&send_req);

	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];

	//Collect sub results
	for (int i = 1; i < n_procs; i++) {
	    double * out_i = &out[n_rows_per_worker*(i-1)];
	    if (i == n_procs-1) {
		MPI_Irecv(out_i, n_rows-n_rows_per_worker*(i-1), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else {
		MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	}

	//Wait for all workers to finish (in code, to measure roundtrip times)
	int n_done = 0, completed[n_workers];
	double worker_stats[n_workers][3];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			worker_stats[i][0] = (double)get_time() - out[n_rows_per_worker*(i+1)-1];
			worker_stats[i][1] = out[n_rows_per_worker*(i+1)-2];
			worker_stats[i][2] = out[n_rows_per_worker*(i+1)-3];
			completed[i] = 1;
			n_done += 1;
		    }
		}
	    }
	}

	int t2 = get_time();

	cout << "-------------------------------------------" << endl;
	cout << "UNCODED MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;
	//Print out worker stats
	for (int i = 1; i < n_procs; i++) {
	    double time_send_solution = worker_stats[i-1][0];
	    double time_receive_input = worker_stats[i-1][1];
	    double work_time = worker_stats[i-1][2];
	    cout << "WORKER " << i << endl << " RECEIVE INPUT TIME: " << time_receive_input << endl << " WORK TIME: " << work_time << endl << " SEND RESULT TIME: " << time_send_solution << endl;
	}
	cout << "-------------------------------------------" << endl;

	free(vec);
	free(out);
    }
    else {
	//Receive the  vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//The last element of the input vector is the time that it was sent
	long long int time_sent = (long long int)subvec[m_cols-1];
	long long int time_received = get_time();
	long long int time_to_send_to_worker = time_received - time_sent;
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	long long int time_to_compute = get_time() - time_received;
	long long int time_resent = get_time();
	//Last element in returned vector is the time to send the solution
	subout[n_rows_per_worker-1] = (double)time_resent;
	//Second to last element in returned vector is the time to receive the input
	subout[n_rows_per_worker-2] = (double)time_to_send_to_worker;
	//Third to last element in returned vector is the time for the worker to do the computation
	subout[n_rows_per_worker-3] = (double)time_to_compute;
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}
