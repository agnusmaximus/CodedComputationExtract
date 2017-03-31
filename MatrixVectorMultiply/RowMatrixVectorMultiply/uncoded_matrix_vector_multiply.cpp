#include <iostream>
#include <mpi.h>
#include "../util.cpp"
#include <cstdlib>
#include <string.h>

void uncoded_matrix_vector_multiply_parallel(int m_cols, int n_rows) {
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
	long long int times[n_workers];

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	    times[i-1] = get_time();
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
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			completed[i] = 1;
			n_done += 1;
			times[i] = get_time() - times[i];
		    }
		}
	    }
	}

	long long int t2 = get_time();

	check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << "UNCODED MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;
	print_worker_roundtrip_times(times, n_workers);
	cout << "-------------------------------------------" << endl;

	free(vec);
	free(out);
    }
    else {

	//Receive the submatrix and vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}

//void uncoded_matrix_vector_multiply_parallel_singlerun(int m_cols, int n_rows) {
void uncoded_matrix_vector_multiply_parallel_singlerun(int m_cols, int n_rows, double *submat, double *subvec, double *subout, double *vec, double *out) {

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

    //Main parallel code
    if (proc_id == 0) {

	//load_input_vector_no_alloc(m_cols, vec);

	//Distribute vectors
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
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
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			completed[i] = 1;
			n_done += 1;
		    }
		}
	    }
	}
    }
    else {

	//Receive the submatrix and vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}


void uncoded_matrix_vector_multiply_parallel_multirun(int m_cols, int n_rows, int num_runs) {
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / n_workers;
    double *submat, *subvec, *subout, *vec, *out;

    //Allocate data
    if (proc_id != 0) {
	submat = (double *)malloc(sizeof(double) * n_rows_per_worker * m_cols);
	load_submatrix_no_alloc(proc_id, m_cols, n_rows_per_worker, submat);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }
    else {
	vec = generate_empty_vector(m_cols);
	out = generate_empty_vector(n_rows);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    long long int t1 = get_time();
    for (int i = 0; i < num_runs; i++) {
	uncoded_matrix_vector_multiply_parallel_singlerun(m_cols, n_rows, submat, subvec, subout, vec, out);
	if (proc_id == 0) {
	    memset(out, 0, sizeof(double) * n_rows);
	}
	else {
	    memset(subout, 0, sizeof(double) * n_rows_per_worker);
	}
    }
    long long int t2 = get_time();
    if (proc_id == 0)
	cout << "ELAPSED TIME: " << t2-t1 << endl;

    //Deallocate data
    if (proc_id != 0) {
	free(submat);
	free(subvec);
	free(subout);
    }
    else {
	free(vec);
	free(out);
    }
}
