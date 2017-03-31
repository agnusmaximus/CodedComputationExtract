#include <iostream>
#include <mpi.h>
#include "../util.cpp"
#include <cstdlib>
#include <string.h>

void uncoded_column_matrix_vector_multiply_parallel(int m_cols, int n_rows) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_cols_per_worker = m_cols / n_workers;

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, n_cols_per_worker, n_rows);
	subvec = (double *)malloc(sizeof(double) * n_cols_per_worker);
	subout = generate_empty_vector(n_rows);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Main parallel code
    if (proc_id == 0) {

	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows);
	double *worker_out = generate_empty_vector(n_rows * n_workers);
	long long int times[n_workers];

	//Distribute vectors
	long long int t1 = get_time();
	MPI_Request send_req;
	for (int i = 1; i < n_procs; i++) {
	    double *partial_vec = &vec[(i-1)*n_cols_per_worker];
	    MPI_Isend(partial_vec, n_cols_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	    times[i-1] = get_time();
	}
	MPI_Request_free(&send_req);

	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];

	//Collect sub results
	for (int i = 1; i < n_procs; i++) {
	    double * out_i = &worker_out[n_rows*(i-1)];
	    MPI_Irecv(out_i, n_rows, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
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
			vector_vector_sum(&worker_out[i*n_rows], out, out, n_rows);
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
	MPI_Recv(subvec, n_cols_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, n_cols_per_worker, n_rows);
	MPI_Send(subout, n_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}
