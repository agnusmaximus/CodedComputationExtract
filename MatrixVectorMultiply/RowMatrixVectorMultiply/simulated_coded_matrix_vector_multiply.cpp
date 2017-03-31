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

void simulated_coded_matrix_vector_multiply_parallel(int m_cols, int n_rows, int parity) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-parity);

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
	double *out = generate_empty_vector(n_rows_per_worker * n_workers);
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
	    MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	}

	//Wait for n-1 workers to complete, computing partial sums as they arrive
	int n_done = 0, completed[n_workers];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_done < n_workers-parity) { //CHANGE BACK
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

	long long int t2 = get_time();
	//check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << parity << " CODED MATRIX VECTOR MULTIPLY\nELAPSED TIME: " << t2-t1 << endl;

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
