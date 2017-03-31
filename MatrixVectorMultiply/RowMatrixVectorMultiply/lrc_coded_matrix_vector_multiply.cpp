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

void decode_missing_split(int n_workers, int n_rows_per_worker,
			      int *completed, double *out,
			      double *rest_sum1, double *rest_sum2, int n_done) {
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
}

// Does 2 coded matrix vector multiply on n_rows/n_splitsxm_cols matrices
void lrc_matrix_vector_multiply_parallel(int m_cols, int n_rows, int n_splits) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_workers_per_split = n_workers  / n_splits;
    int n_rows_per_split = n_rows / n_splits;
    int n_rows_per_worker = (n_rows / n_splits) / (n_workers_per_split-2);

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout, *parity1, *parity2;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (proc_id == 0) {
	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows);
	double *rest_sum1[n_splits], *rest_sum2[n_splits];
	double *parity1[n_splits], *parity2[n_splits];
	for (int i = 0; i < n_splits; i++) {
	    rest_sum1[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	    rest_sum2[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	    parity1[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	    parity2[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	}

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
	    int n_parity_skipped = (i-1)/n_workers_per_split * 2;
	    double * out_i = &out[n_rows_per_worker*(i-1-n_parity_skipped)];
	    if ((i-1) % n_workers_per_split == n_workers_per_split-2) {
		MPI_Irecv(parity1[(i-1)/n_workers_per_split], n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else if ((i-1) % n_workers_per_split == n_workers_per_split-1) {
		MPI_Irecv(parity2[(i-1)/n_workers_per_split], n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else {
		MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	}

	//Wait for all splits to finish
	int splits_completed[n_splits], n_splits_completed = 0;
	int n_workers_completed_per_split[n_splits];
	int completed[n_workers], n_workers_done = 0;
	for (int i = 0; i < n_splits; i++) n_workers_completed_per_split[i] = splits_completed[i] = 0;
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_splits_completed != n_splits) {
	    for (int i = 0; i < n_workers; i++) {
		int cur_split = i/n_workers_per_split;
		int n_parity_skipped = i/n_workers_per_split * 2;
		//Only check worker if not completed the split the worker was working on
		if (splits_completed[cur_split] == 0 && !completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			n_workers_done++;
			times[i] = get_time() - times[i];
			completed[i] = 1;

			//Receive partial results
			if (i % n_workers_per_split == n_workers_per_split - 2) {
			    vector_vector_sum(rest_sum1[cur_split], parity1[cur_split], rest_sum1[cur_split], n_rows_per_worker);
			}
			else if (i % n_workers_per_split == n_workers_per_split - 1) {
			    vector_vector_sum(rest_sum2[cur_split], parity2[cur_split], rest_sum2[cur_split], n_rows_per_worker);
			}
			else {
			    vector_vector_subtract(rest_sum1[cur_split], &out[n_rows_per_worker*(i-n_parity_skipped)], rest_sum1[cur_split], n_rows_per_worker);
			    vector_vector_add_scalark(rest_sum2[cur_split], &out[n_rows_per_worker*(i-n_parity_skipped)], rest_sum2[cur_split], -(i%n_workers_per_split+1), n_rows_per_worker);
			}

			n_workers_completed_per_split[cur_split]++;
		    }
		}

		//Another pass to potentially decode individual splits
		if (splits_completed[cur_split] == 0 && n_workers_completed_per_split[cur_split] >= n_workers_per_split-2) {
		    decode_missing_split(n_workers_per_split, n_rows_per_worker,
					 &completed[cur_split*n_workers_per_split], &out[cur_split*n_rows_per_split],
					 rest_sum1[cur_split], rest_sum2[cur_split],
					 n_workers_completed_per_split[cur_split]);
		    n_splits_completed++;
		    splits_completed[cur_split] = 1;
		}
	    }
	}

	long long int t2 = get_time();
	check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << "LRC 2  " << n_splits << " SPLIT MATRIX VECTOR MULTIPLY\nELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	while (n_workers_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_workers_done += 1;
		    }
		}
	    }
	}
	print_worker_roundtrip_times(times, n_workers);
	cout << "-------------------------------------------" << endl;

	for (int i = 0; i < n_splits; i++) {
	    free(rest_sum1[i]);
	    free(rest_sum2[i]);
	    free(parity1[i]);
	    free(parity2[i]);
	}
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

void lrc_coded_matrix_vector_multiply_parallel(int m_cols, int n_rows, int n_splits) {

    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_workers_per_split = n_workers  / n_splits;
    int n_rows_per_split = n_rows / (n_splits-2);
    int n_rows_per_worker = n_rows_per_split / (n_workers_per_split-2);

    //Preprocessing - generate random sub-matrices (for simulating
    //the matrix already pre-loaded on the worker)
    double *submat, *subvec, *subout, *parity1, *parity2;
    if (proc_id != 0) {
	submat = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	subvec = (double *)malloc(sizeof(double) * m_cols);
	subout = generate_empty_vector(n_rows_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (proc_id == 0) {
	double *vec = load_input_vector(m_cols);
	double *out = generate_empty_vector(n_rows+2*n_rows_per_split);
	double *rest_sum1[n_splits], *rest_sum2[n_splits];
	double *parity1[n_splits], *parity2[n_splits];
	for (int i = 0; i < n_splits; i++) {
	    rest_sum1[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	    rest_sum2[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	    parity1[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	    parity2[i] = (double *)calloc(n_rows_per_worker, sizeof(double));
	}
	double *rest_split_sum1 = (double *)calloc(n_rows_per_split, sizeof(double));
	double *rest_split_sum2 = (double *)calloc(n_rows_per_split, sizeof(double));

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
	    int n_parity_skipped = (i-1)/n_workers_per_split * 2;
	    double * out_i = &out[n_rows_per_worker*(i-1-n_parity_skipped)];
	    int cur_split = (i-1)/n_workers_per_split;
	    if ((i-1) % n_workers_per_split == n_workers_per_split-2) {
		MPI_Irecv(parity1[cur_split], n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else if ((i-1) % n_workers_per_split == n_workers_per_split-1) {
		MPI_Irecv(parity2[cur_split], n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	    else {
		MPI_Irecv(out_i, n_rows_per_worker, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_reqs[i-1]);
	    }
	}

	//Wait for all splits to finish
	int splits_completed[n_splits], n_splits_completed = 0;
	int n_workers_completed_per_split[n_splits];
	int completed[n_workers], n_workers_done = 0;
	for (int i = 0; i < n_splits; i++) n_workers_completed_per_split[i] = splits_completed[i] = 0;
	for (int i = 0; i < n_workers; i++) completed[i] = 0;
	while (n_splits_completed < n_splits-2) {
	    for (int i = 0; i < n_workers; i++) {
		int cur_split = i/n_workers_per_split;
		int n_parity_skipped = i/n_workers_per_split * 2;
		//Only check worker if not completed the split the worker was working on
		if (splits_completed[cur_split] == 0 && !completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			n_workers_done++;
			times[i] = get_time() - times[i];
			completed[i] = 1;

			//Receive partial results
			if (i % n_workers_per_split == n_workers_per_split - 2) {
			    vector_vector_sum(rest_sum1[cur_split], parity1[cur_split], rest_sum1[cur_split], n_rows_per_worker);
			}
			else if (i % n_workers_per_split == n_workers_per_split - 1) {
			    vector_vector_sum(rest_sum2[cur_split], parity2[cur_split], rest_sum2[cur_split], n_rows_per_worker);
			}
			else {
			    vector_vector_subtract(rest_sum1[cur_split], &out[n_rows_per_worker*(i-n_parity_skipped)], rest_sum1[cur_split], n_rows_per_worker);
			    vector_vector_add_scalark(rest_sum2[cur_split], &out[n_rows_per_worker*(i-n_parity_skipped)], rest_sum2[cur_split], -(i%n_workers_per_split+1), n_rows_per_worker);
			}

			n_workers_completed_per_split[cur_split]++;
		    }
		}

		//Another pass to potentially decode individual splits
		if (splits_completed[cur_split] == 0 && n_workers_completed_per_split[cur_split] >= n_workers_per_split-2) {
		    decode_missing_split(n_workers_per_split, n_rows_per_worker,
					 &completed[cur_split*n_workers_per_split], &out[cur_split*n_rows_per_split],
					 rest_sum1[cur_split], rest_sum2[cur_split],
					 n_workers_completed_per_split[cur_split]);
		    n_splits_completed++;
		    splits_completed[cur_split] = 1;



		    //Parity splits
		    if (cur_split == n_splits - 2) {
			vector_vector_sum(rest_split_sum1, &out[n_rows_per_split*cur_split], rest_split_sum1, n_rows_per_split);
		    }
		    else if (cur_split == n_splits - 1) {
			vector_vector_sum(rest_split_sum2, &out[n_rows_per_split*cur_split], rest_split_sum2, n_rows_per_split);
		    }
		    else {
			vector_vector_subtract(rest_split_sum1, &out[n_rows_per_split*(cur_split)], rest_split_sum1, n_rows_per_split);
			vector_vector_add_scalark(rest_split_sum2, &out[n_rows_per_split*(cur_split)], rest_split_sum2, -(cur_split+1), n_rows_per_split);
		    }
		}
	    }
	}

	//Decode last splits
	decode_missing_split(n_splits, n_rows_per_split, splits_completed, out,
			     rest_split_sum1, rest_split_sum2,
			     n_splits_completed);

	long long int t2 = get_time();
	check_correct(out, n_rows);

	cout << "-------------------------------------------" << endl;
	cout << "TWO CODED " << n_splits << " SPLIT MATRIX VECTOR MULTIPLY\nELAPSED TIME: " << t2-t1 << endl;

	//Wait for last process to finish
	while (n_workers_done < n_workers) {
	    for (int i = 0; i < n_workers; i++) {
		if (!completed[i]) {
		    int did_complete = 0;
		    MPI_Test(&recv_reqs[i], &did_complete, MPI_STATUS_IGNORE);
		    if (did_complete) {
			times[i] = get_time() - times[i];
			completed[i] = 1;
			n_workers_done += 1;
		    }
		}
	    }
	}
	print_worker_roundtrip_times(times, n_workers);
	cout << "-------------------------------------------" << endl;

	for (int i = 0; i < n_splits; i++) {
	    free(rest_sum1[i]);
	    free(rest_sum2[i]);
	    free(parity1[i]);
	    free(parity2[i]);
	}
	free(vec);
	free(out);
	free(rest_split_sum1);
	free(rest_split_sum2);

    }
    else {
	//if (proc_id >= 5 && proc_id <= 12) sleep(5);
	//Receive the  vector
	MPI_Recv(subvec, m_cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	matrix_vector_multiply(submat, subvec, subout, m_cols, n_rows_per_worker);
	MPI_Send(subout, n_rows_per_worker, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	free(submat);
	free(subvec);
	free(subout);
    }
}
