#include <iostream>
#include <mpi.h>
#include "../MatrixVectorMultiply/util.cpp"
#include "LinearRegressionRoutines.cpp"
#include "LinearRegressionRoutinesV2.cpp"
#include <cstdlib>
#include <math.h>

using namespace std;

#define UNCODED_LINEAR_REGRESSION_VAL 0
#define CODED_LINEAR_REGRESSION_VAL 1

double calculate_error(double *out, int n_procs) {
    int n_workers = n_procs - 1;
    int n_rows_per_worker = 0;
    if (TYPE == UNCODED_LINEAR_REGRESSION_VAL) {
	n_rows_per_worker = N_ROWS / n_workers;
    }
    else {
	n_rows_per_worker = N_ROWS / (n_workers-2);
    }
    double *mat = (double *)malloc(sizeof(double) * N_ROWS * M_COLS);
    int n_rows_copied_so_far = 0;
    for (int i = 1; i < n_procs; i++) {
	if (n_rows_copied_so_far >= N_ROWS) break;
	load_submatrix_no_alloc(i, M_COLS, n_rows_per_worker, &mat[(i-1) * (M_COLS * n_rows_per_worker)]);
	n_rows_copied_so_far += n_rows_per_worker;
    }
    double *y = load_vector(N_ROWS, "labels.dat");
    double *ax = (double *)calloc(N_ROWS, sizeof(double));
    matrix_vector_multiply(mat, out, ax, M_COLS, N_ROWS);
    double ax_min_y = 0;
    double y2 = 0;
    for (int i = 0; i < N_ROWS; i++) {
	ax_min_y += (ax[i]-y[i]) * (ax[i]-y[i]);
	y2 += y[i] * y[i];
    }
    free(mat);
    free(y);
    free(ax);
    return (ax_min_y / y2) / N_ROWS;
}

int main(void) {
    //Initialization
    //Initialization
    srand(time(NULL));
    int n_procs, proc_id;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Barrier(MPI_COMM_WORLD);

    double *out;
    long long int t1, t2;
    if (TYPE == UNCODED_LINEAR_REGRESSION_VAL) {
	t1 = get_time();
	out = uncoded_linear_regression(M_COLS, N_ROWS, DATA_PATH, LEARNING_RATE);
	t2 = get_time();
    }
    else {
	t1 = get_time();
	out = coded2_linear_regression_v2(M_COLS, N_ROWS, DATA_PATH, LEARNING_RATE);
	t2 = get_time();
    }
    if (proc_id == 0) {
	cout << "ELAPSED TIME: " << t2-t1 << endl;
	for (int i = 0; i < min(10, M_COLS); i++)
	    cout << out[i] << " ";
	cout << endl;

	//Compute msq
	double msq = calculate_error(out, n_procs);
	cout << "ERROR: " << msq << endl;
    }
    //cout << "DONE " << proc_id << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_id == 0)
      free(out);
    MPI_Finalize();
}
