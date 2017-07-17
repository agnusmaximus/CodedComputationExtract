#include <iostream>
#include <mpi.h>
#include "../MatrixVectorMultiply/RowMatrixVectorMultiply/uncoded_matrix_vector_multiply.cpp"
#include "../MatrixVectorMultiply/RowMatrixVectorMultiply/coded_matrix_vector_multiply.cpp"
#include "../MatrixVectorMultiply/util.cpp"

using namespace std;

#define THRESHOLD .001

double error(double *diffs, int n) {
    double val = 0;
    for (int i = 0; i < n; i++) {
	val += diffs[i] * diffs[i];
    }
    return val;
}

bool is_done(double *gradient, double lr, int n) {
    double diff = 0;
    for (int i = 0; i < n; i++) {
	double val = gradient[i] * lr;
	diff += val * val;
    }
    return diff < THRESHOLD;
}

double * coded2_linear_regression(int m_cols, int n_rows, string data_path, double learning_rate) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / (n_workers-2);
    int m_cols_per_worker = m_cols / (n_workers-2);

    //Master process data
    double *out_transpose, *gradient, *out, *outkk, *outkkt;
    double *labels;
    double *rest_sum1, *rest_sum2, *parity1, *parity2;

    //Worker process data
    double *sub_data;
    double *sub_vec, *sub_out;
    double *sub_data_transpose;
    double *sub_vec_transpose, *sub_out_transpose;

    int max_row_col = max(n_rows_per_worker, m_cols_per_worker);

    //Load memory
    if (proc_id == 0) {
	gradient = (double *)calloc(m_cols, sizeof(double));
	out_transpose = (double *)calloc(n_rows, sizeof(double));
	out = (double *)calloc(m_cols, sizeof(double));
	labels = load_vector(n_rows, "labels.dat");

	rest_sum1 = (double *)calloc(max_row_col, sizeof(double));
	rest_sum2 = (double *)calloc(max_row_col, sizeof(double));
	parity1 = (double *)calloc(max_row_col, sizeof(double));
	parity2 = (double *)calloc(max_row_col, sizeof(double));
	outkk = (double *)calloc(max(n_rows, m_cols), sizeof(double));
	outkkt = (double *)calloc(max(n_rows, m_cols), sizeof(double));
    }
    else {
	sub_vec = (double *)calloc(m_cols, sizeof(double));
	sub_out = (double *)calloc(n_rows_per_worker, sizeof(double));
	sub_vec_transpose = (double *)calloc(n_rows, sizeof(double));
	sub_out_transpose = (double *)calloc(m_cols_per_worker, sizeof(double));
	sub_data = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	sub_data_transpose = load_submatrix_t(proc_id, n_rows, m_cols_per_worker);
    }

    //Do a single gradient descent iteration before the loop
    for (int i = 0; i < N_ITERS; i++) {

	//out = x^(t), so out_transpose = Ax^(t)
	//uncoded_matrix_vector_multiply_parallel_singlerun(m_cols, n_rows, sub_data, sub_vec, sub_out, out, out_transpose);
	coded2_matrix_vector_multiply_parallel_singlerun(m_cols, n_rows, i,
							 sub_data, sub_vec, sub_out,
							 out, out_transpose, outkk,
							 rest_sum1, rest_sum2,
							 parity1, parity2);

	if (proc_id == 0) {
	    //out_transpose = out_transpose - y = Ax^(t) - y
	    vector_vector_subtract(out_transpose, labels, out_transpose, n_rows);

	    memset(rest_sum1, 0, sizeof(double) * max_row_col);
	    memset(rest_sum2, 0, sizeof(double) * max_row_col);
	    memset(parity1, 0, sizeof(double) * max_row_col);
	    memset(parity2, 0, sizeof(double) * max_row_col);
	}

	//gradient = A^T(Ax^(t) - y)
	//uncoded_matrix_vector_multiply_parallel_singlerun(n_rows, m_cols, sub_data_transpose, sub_vec_transpose, sub_out_transpose, out_transpose, gradient);
	coded2_matrix_vector_multiply_parallel_singlerun(n_rows, m_cols, i + N_ITERS,
							 sub_data_transpose, sub_vec_transpose, sub_out_transpose,
							 out_transpose, gradient, outkkt,
							 rest_sum1, rest_sum2,
							 parity1, parity2);
	//MPI_Barrier(MPI_COMM_WORLD);
	//break;


       if (proc_id == 0) {
	   //out = out - learning_rate * gradient
	   vector_vector_subtract_scalark(out, gradient, out, learning_rate, m_cols);
       }

	if (proc_id == 0) {
	    memset(out_transpose, 0, sizeof(double) * n_rows);
	    memset(gradient, 0, sizeof(double) * m_cols);

	    memset(rest_sum1, 0, sizeof(double) * max_row_col);
	    memset(rest_sum2, 0, sizeof(double) * max_row_col);
	    memset(parity1, 0, sizeof(double) * max_row_col);
	    memset(parity2, 0, sizeof(double) * max_row_col);
	}
	else {
	    memset(sub_out, 0, sizeof(double) * n_rows_per_worker);
	    memset(sub_out_transpose, 0, sizeof(double) * m_cols_per_worker);
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (proc_id == 0) {
	free(outkk);
	free(gradient);
	free(out_transpose);
	free(labels);
	free(rest_sum1);
	free(rest_sum2);
	free(parity1);
	free(parity2);
	free(outkkt);
    }
    else {
	free(sub_data);
	free(sub_vec);
	free(sub_out);
	free(sub_data_transpose);
	free(sub_vec_transpose);
	free(sub_out_transpose);
    }
    return out;
}

double * uncoded_linear_regression(int m_cols, int n_rows, string data_path, double learning_rate) {
    //Initialization
    int n_procs, proc_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Get_processor_name(processor_name, &name_len);
    int n_workers = n_procs-1;
    int n_rows_per_worker = n_rows / n_workers;
    int m_cols_per_worker = m_cols / n_workers;

    //Master process data
    double *out_transpose, *gradient, *out;
    double *labels;

    //Worker process data
    double *sub_data;
    double *sub_vec, *sub_out;
    double *sub_data_transpose;
    double *sub_vec_transpose, *sub_out_transpose;

    //Load memory
    if (proc_id == 0) {
	gradient = (double *)calloc(m_cols, sizeof(double));
	out_transpose = (double *)calloc(n_rows, sizeof(double));
	out = (double *)calloc(m_cols, sizeof(double));
	labels = load_vector(n_rows, "labels.dat");
    }
    else {
	sub_vec = (double *)calloc(m_cols, sizeof(double));
	sub_out = (double *)calloc(n_rows_per_worker, sizeof(double));
	sub_vec_transpose = (double *)calloc(n_rows, sizeof(double));
	sub_out_transpose = (double *)calloc(m_cols_per_worker, sizeof(double));
	sub_data = load_submatrix(proc_id, m_cols, n_rows_per_worker);
	sub_data_transpose = load_submatrix_t(proc_id, n_rows, m_cols_per_worker);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Do a single gradient descent iteration before the loop
    for (int i = 0; i < N_ITERS; i++) {

	//out = x^(t), so out_transpose = Ax^(t)
	uncoded_matrix_vector_multiply_parallel_singlerun(m_cols, n_rows, sub_data, sub_vec, sub_out, out, out_transpose);

	if (proc_id == 0) {
	    //out_transpose = out_transpose - y = Ax^(t) - y
	    vector_vector_subtract(out_transpose, labels, out_transpose, n_rows);
	}
	//gradient = A^T(Ax^(t) - y)
	uncoded_matrix_vector_multiply_parallel_singlerun(n_rows, m_cols, sub_data_transpose, sub_vec_transpose, sub_out_transpose, out_transpose, gradient);

	if (proc_id == 0) {
	    //out = out - learning_rate * gradient
	    vector_vector_subtract_scalark(out, gradient, out, learning_rate, m_cols);
	}

	if (proc_id == 0) {
	    memset(out_transpose, 0, sizeof(double) * n_rows);
	    memset(gradient, 0, sizeof(double) * m_cols);
	}
	else {
	    memset(sub_out, 0, sizeof(double) * n_rows_per_worker);
	    memset(sub_out_transpose, 0, sizeof(double) * m_cols_per_worker);
	}
    }

    if (proc_id == 0) {
	free(gradient);
	free(out_transpose);
	free(labels);
    }
    else {
	free(sub_data);
	free(sub_vec);
	free(sub_out);
	free(sub_data_transpose);
	free(sub_vec_transpose);
	free(sub_out_transpose);
    }
    return out;
}
