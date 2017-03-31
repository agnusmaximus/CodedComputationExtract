#include <iostream>
#include "../util.cpp"
#include <mpi.h>
#include "uncoded_column_matrix_vector_multiply.cpp"

#define UNCODED_VAL 0

int main(void) {
    //Initialization
    srand(time(NULL));
    int n_procs, proc_id;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    MPI_Barrier(MPI_COMM_WORLD);

    if (CODED == UNCODED_VAL) {
	uncoded_column_matrix_vector_multiply_parallel(N_COLS, N_ROWS);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
