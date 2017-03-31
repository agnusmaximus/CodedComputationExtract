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

void replicated_matrix_vector_multiply_parallel_profile(int m_cols, int n_rows, int n_chunks, int n_workers_per_chunk) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / n_chunks;

    //Preprocessing - read random sub-matrices (for simulating
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
	double *sub_out = generate_empty_vector(n_rows * n_workers);
	double *out = generate_empty_vector(n_rows);
	long long int times[n_workers];

	//Distribute input vectors to all workers
	MPI_Request send_req;
	long long int t1 = get_time();
	for (int i = 1; i < n_procs; i++) {
	    times[i-1] = get_time();
	    MPI_Isend(vec, m_cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &send_req);
	}
	MPI_Request_free(&send_req);

	//Collect results
	MPI_Request recv_reqs[n_workers];
	MPI_Status statuses[n_workers];
	for (int i = 1; i < n_procs; i++) {
	    double * out_i = &sub_out[n_rows_per_worker * (i-1)];
	    MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	}

	int distinct_chunks_collected[n_chunks];
	int num_distinct_chunks_collected = 0, n_done = 0;
	for (int i = 0; i < n_chunks; i++) distinct_chunks_collected[i] = 0;
	int completed[n_workers];
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (num_distinct_chunks_collected != n_chunks) {
	    //Check statuses of all workers
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);

		    //This worker has completed is sub matrix multiplication
		    if (did_complete) {
			completed[i] = 1;
			n_done += 1;
			times[i] = get_time() - times[i];

			//Do we already have his result? (another worker doing same chunk already finished?)
			//If not, copy over to final result
			if (distinct_chunks_collected[i/n_workers_per_chunk] == 0) {
			    distinct_chunks_collected[i/n_workers_per_chunk] = 1;
			    double *sub_out_i = &sub_out[n_rows_per_worker*i];
			    double *out_i = &out[i/n_workers_per_chunk * n_rows_per_worker];
			    memcpy(out_i, sub_out_i, sizeof(double) * n_rows_per_worker);
			    num_distinct_chunks_collected++;
			}
		    }
		}
	    }
	}
	long long int t2 = get_time();
	check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << "REPLICATED CODED MATRIX VECTOR MULTIPLY (" << n_chunks << " CHUNKS, " << n_workers_per_chunk << " WORKERS PER CHUNK)\n ELAPSED TIME: " << t2-t1 << endl;

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
	free(sub_out);
	free(out);
    }
    else {
	//Worker stuff
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}
