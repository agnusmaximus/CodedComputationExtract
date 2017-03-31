#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <cblas.h>
#include "../util.cpp"
#include "uncoded_matrix_vector_multiply.cpp"
#include "coded_matrix_vector_multiply.cpp"
#include "replicated_matrix_vector_multiply.cpp"
#include "coded_matrix_vector_multiply_profiling.cpp"
#include "simulated_coded_matrix_vector_multiply.cpp"
#include "lrc_coded_matrix_vector_multiply.cpp"

#ifndef PROFILE
#define PROFILE 0
#endif

#define SERIAL_BENCH 0
#define LRC_CODED_VAL -5
#define LRC_VAL -4
#define SIMULATED_CODED_VAL -3
#define REPLICATION_VAL -2
#define EXTREME_VAL -1
#define UNCODED_VAL 0
#define CODED_VAL 1
#define CODED2_VAL 2
#define CODED2RUN_VAL 3
#define UNCODEDRUN_VAL 4

using namespace std;

void serial_benchmark(int m_cols, int n_rows) {
    int proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    if (proc_id == 0) {
	double *mat = generate_random_matrix(m_cols, n_rows);
	double *vec = generate_random_vector(m_cols);
	double *out = generate_empty_vector(n_rows);

	long long int t1 = get_time();
	matrix_vector_multiply_naive(mat, vec, out, m_cols, n_rows);
	long long int t2 = get_time();
	cout << "SERIAL MATRIX VECTOR MULTIPLY ELAPSED TIME: " << t2-t1 << endl;
	free(mat);
	free(vec);
	free(out);
    }
}

int main(void) {
    //Initialization
    srand(time(NULL));
    int n_procs, proc_id;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    MPI_Barrier(MPI_COMM_WORLD);

    if (PROFILE) {
	if (CODED == CODED_VAL) {
	    coded_matrix_vector_multiply_parallel_profile(N_COLS, N_ROWS);
	}
	if (CODED == CODED2_VAL) {
	    coded2_matrix_vector_multiply_parallel_profile(N_COLS, N_ROWS);
	}
	if (CODED == UNCODED_VAL) {
	    uncoded_matrix_vector_multiply_parallel_profile(N_COLS, N_ROWS);
	}
    }
    else {
	if (CODED == LRC_VAL) {
	    #ifdef N_SPLITS
	    lrc_matrix_vector_multiply_parallel(N_COLS, N_ROWS, N_SPLITS);
	    #endif
	}
	if (CODED == LRC_CODED_VAL) {
	    #ifdef N_SPLITS
	    lrc_coded_matrix_vector_multiply_parallel(N_COLS, N_ROWS, N_SPLITS);
	    #endif
	}
	if (CODED == REPLICATION_VAL) {
            #if defined(N_WORKERS_PER_CHUNK) && defined(N_CHUNKS)
	    replicated_matrix_vector_multiply_parallel_profile(N_COLS, N_ROWS, N_CHUNKS, N_WORKERS_PER_CHUNK);
            #endif
	}
	if (CODED == SIMULATED_CODED_VAL) {
	    #ifdef PARITY
	    simulated_coded_matrix_vector_multiply_parallel(N_COLS, N_ROWS, PARITY);
	    #endif
	}
	if (CODED == UNCODED_VAL) {
	    uncoded_matrix_vector_multiply_parallel(N_COLS, N_ROWS);
	}
	else if (CODED == CODED_VAL) {
	    coded_matrix_vector_multiply_parallel(N_COLS, N_ROWS);
	}
	else if (CODED == CODED2_VAL) {
	    coded2_matrix_vector_multiply_parallel(N_COLS, N_ROWS);
	}
	else if (CODED == EXTREME_VAL) {
	    extreme_matrix_vector_multiply_parallel(N_COLS, N_ROWS);
	}
	else if (CODED == CODED2RUN_VAL) {
	    coded2_matrix_vector_multiply_parallel_multirun(N_COLS, N_ROWS, NUM_RUNS);
	}
	else if (CODED == UNCODEDRUN_VAL) {
	    uncoded_matrix_vector_multiply_parallel_multirun(N_COLS, N_ROWS, NUM_RUNS);
	}
	if (SERIAL_BENCH) serial_benchmark(N_COLS, N_ROWS);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
