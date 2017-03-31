#ifndef UTIL
#define UTIL

#include <time.h>
#include <iostream>
#include <cblas.h>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

template <typename T>
string int_to_string(T Number) {
    ostringstream ss;
    ss << Number;
    return ss.str();
}

double rand_double() {
    return rand() % 100;
}

long long int get_time() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

void verify_equivalent(double *a, double *b, int len) {
    for (int i = 0; i < len; i++)
	if (a[i] != b[i]) {
	    cout << "ERROR!" << endl;
	    return;
	}
    cout << "ALL CORRECT!" << endl;
}

void print_matrix(double *mat, int m_cols, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	cout << " | ";
	for (int j = 0; j < m_cols; j++) {
	    cout << mat[i*m_cols + j] << " ";
	}
	cout << " | " << endl;
    }
    cout << endl;
}

void print_vect(double *vec, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	cout << " | " << vec[i] << " | " << endl;
    }
    cout << endl;
}

double * generate_random_matrix(int m_cols, int n_rows) {
    double * mat = (double *)malloc(sizeof(double) * m_cols * n_rows);
    for (int i = 0; i < m_cols * n_rows; i++)
	mat[i] = rand_double();
    return mat;
}

double * generate_random_vector(int n_rows) {
    double * vect = (double *)malloc(sizeof(double) * n_rows);
    for (int i = 0; i < n_rows; i++)
	vect[i] = rand_double();
    return vect;
}

double * generate_empty_vector(int n_rows) {
    double * vect = (double *)calloc(n_rows, sizeof(double));
    return vect;
}

void matrix_vector_multiply(double *mat, double *vec, double *out, int m_cols, int n_rows) {
    cblas_dgemv(CblasRowMajor,  CblasNoTrans,  n_rows, m_cols, 1, mat, m_cols, vec, 1, 1, out, 1);
}

void matrix_matrix_sum(double *mat1, double *mat2, double *out, int m_cols, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	for (int j = 0; j < m_cols; j++) {
	    out[i*m_cols+j] = mat1[i*m_cols+j] + mat2[i*m_cols+j];
	}
    }
}

void vector_vector_scalar_multiply(double *v1, double *out, double k, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	out[i] = v1[i] * k;
    }
}

void vector_vector_sum(double *v1, double *v2, double *out, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	out[i] = v1[i] + v2[i];
    }
}

void vector_vector_subtract(double *v1, double *v2, double *out, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	out[i] = v1[i] - v2[i];
    }
}

void vector_vector_subtract_scalark(double *v1, double *v2, double *out, double k, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	out[i] = v1[i] - v2[i]*k;
    }
}

void vector_vector_add_scalark(double *v1, double *v2, double *out, double k, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	out[i] = v1[i] + v2[i]*k;
    }
}

void matrix_vector_multiply_naive(double *mat, double *vec, double *out, int m_cols, int n_rows) {
    for (int i = 0; i < n_rows; i++) {
	for (int j = 0; j < m_cols; j++) {
	    out[i] += mat[i*m_cols+j] * vec[j];
	}
    }
}

void check(double *mat, double *vec, double *out, int m_cols, int n_rows) {
    double *real_out = generate_empty_vector(n_rows);
    matrix_vector_multiply(mat, vec, real_out, m_cols, n_rows);
    verify_equivalent(real_out, out, n_rows);
}

double * load_submatrix_t(int proc_id, int m_cols, int n_rows) {
    string fname = DATA_PATH + int_to_string(proc_id) + "t.dat";
    ifstream in(fname.c_str());
    if (!in) {
	cout << "ERROR OPENING FILE: " << fname << endl;
	exit(0);
    }
    double *ret = (double *)malloc(sizeof(double) * m_cols * n_rows);
    for (int i = 0; i < n_rows; i++) {
	for (int j = 0; j < m_cols; j++) {
	    int value;
	    in >> value;
	    ret[i*m_cols+j] = value;
	}
    }

    in.close();
    return ret;
}

double  * load_submatrix(int proc_id, int m_cols, int n_rows) {
    string fname = DATA_PATH + int_to_string(proc_id) + ".dat";
    ifstream in(fname.c_str());
    if (!in) {
	cout << "ERROR OPENING FILE: " << fname << endl;
	exit(0);
    }
    double *ret = (double *)malloc(sizeof(double) * m_cols * n_rows);
    for (int i = 0; i < n_rows; i++) {
	for (int j = 0; j < m_cols; j++) {
	    double value;
	    in >> value;
	    ret[i*m_cols+j] = value;
	}
    }

    in.close();
    return ret;
}

void load_submatrix_no_alloc(int proc_id, int m_cols, int n_rows, double *out) {
    string fname = DATA_PATH + int_to_string(proc_id) + ".dat";
    ifstream in(fname.c_str());
    if (!in) {
	cout << "ERROR OPENING FILE: " << fname << endl;

	exit(0);
    }
    for (int i = 0; i < n_rows; i++) {
	for (int j = 0; j < m_cols; j++) {
	    in >> out[i*m_cols+j];
	}
    }

    in.close();
}

double * load_vector(int m, string name) {
    string fname = DATA_PATH + name;
    ifstream in(fname.c_str());
    if (!in) {
	cout << "ERROR OPENING FILE: " << fname << endl;
	exit(0);
    }
    double *ret = (double *)malloc(sizeof(double) * m);
    for (int i = 0; i < m; i++) in >> ret[i];
    in.close();
    return ret;
}

double * load_input_vector(int m) {
    return load_vector(m, "vec.dat");
}

void print_worker_roundtrip_times(long long int *times, int n) {
    for (int i = 0; i < n; i++) {
	cout << "WORKER " << i+1 << " ROUNDTRIP TIME: " << times[i] << endl;
    }
}

void check_correct(double *out, int n) {
    double *correct = load_vector(n, "result.dat");
    for (int i = 0; i < n; i++) {
	if (fabs(correct[i] - out[i]) > 5.96e-08) {
	    cout << "WRONG! " << correct[i] << " vs " << out[i] << ", at index " << i << endl;
	}
    }
}

void load_input_vector_no_alloc(int n, double *out) {
    string fname = DATA_PATH + (string)"vec.dat";
    ifstream in(fname.c_str());
    if (!in) {
	cout << "ERROR OPENING FILE: " << fname << endl;
	exit(0);
    }
    for (int i = 0; i < n; i++) in >> out[i];
    in.close();
}

#endif
