# TODO: Finish the rest

from __future__ import print_function
import sys
import random
import os
from sklearn import linear_model

random.seed(0)

def matrix_vector_multiply(matrix, vector, n_rows, n_cols):
    result = [0 for x in range(n_rows)]
    for i in range(n_rows):
        for j in range(n_cols):
            result[i] += matrix[i][j] * vector[j]
    return result

def generate_empty_matrix(n_rows, n_cols):
    return [[0 for i in range(n_cols)] for x in range(n_rows)]

def scal_mult_matrix_ret(m, s, n_rows, n_cols):
    ret = [[0 for i in range(n_cols)] for j in range(n_rows)]
    for i in range(n_rows):
        for j in range(n_cols):
            ret[i][j] = m[i][j] * s
    return ret

def sum_matrices(m1, m2, n_rows, n_cols):
    # m1 = m1 + m2
    for i in range(n_rows):
        for j in range(n_cols):
            m1[i][j] += m2[i][j]

def save_matrix(m, output):
    f = open(output, "w")
    for i in range(len(m)):
        print(" ".join([str(x) for x in m[i]]), file=f)
    f.close()

def save_vector(m, output):
    f = open(output, "w")
    for i in range(len(m)):
        print(str(m[i]) + " ", file=f)
    f.close()

def generate_random_matrix(n_rows, n_cols):
    return [[random.randint(1, 100) - 50  for i in range(n_cols)] for x in range(n_rows)]

def generate_random_vector(n):
    return [random.randint(1, 100) - 50 for x in range(n)]

if len(sys.argv) != 7:
    print("Usage: python generate_linear_regression_data.py n_procs is_coded n_cols n_rows output_dir should_output_result")
    sys.exit(0)

n_procs, is_coded, n_cols, n_rows, output_dir, should_output_result = sys.argv[1:]
n_procs, is_coded, n_cols, n_rows, should_output_result = (int(x) for x in [n_procs,is_coded,n_cols,n_rows,should_output_result])
output_dir = output_dir + "/" if output_dir[-1] != "/" else output_dir
n_workers, n_rows_per_worker, n_cols_per_worker = -1, -1, -1

n_workers = n_procs-1

if not is_coded:
    n_rows_per_worker = n_rows / n_workers
    n_cols_per_worker = n_cols / n_workers
elif is_coded == 2:
    n_rows_per_worker = n_rows / (n_workers-2)
    n_cols_per_worker = n_cols / (n_workers-2)

if not is_coded:
    assert(n_cols % n_workers == 0)
    assert(n_rows % n_workers == 0)
elif is_coded == 2:
    assert(n_rows % (n_workers-2) == 0)
    assert(n_cols % (n_workers-2) == 0)

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

feature_matrix = generate_random_matrix(n_rows, n_cols)
transpose_matrix = [[feature_matrix[i][j] for i in range(n_rows)] for j in range(n_cols)]
answer_vector = generate_random_vector(n_cols)
labels_vector = matrix_vector_multiply(feature_matrix, answer_vector, n_rows, n_cols)

if not is_coded:
    for i in range(1, n_procs):
        sub_matrix = feature_matrix[(i-1)*n_rows_per_worker:(i)*n_rows_per_worker]
        save_matrix(sub_matrix, output_dir + str(i) + ".dat")
        sub_matrix_transpose = transpose_matrix[(i-1)*n_cols_per_worker:(i)*n_cols_per_worker]
        save_matrix(sub_matrix_transpose, output_dir + str(i) + "t.dat")
    save_vector(labels_vector, output_dir + "labels.dat")
    regr = linear_model.LinearRegression()
    regr.fit(feature_matrix, labels_vector)
    answer = regr.coef_
    save_vector(answer, output_dir + "answer.dat")
if is_coded == 2:
    parity1 = generate_empty_matrix(n_rows_per_worker, n_cols)
    parity2 = generate_empty_matrix(n_rows_per_worker, n_cols)
    parity1t = generate_empty_matrix(n_cols_per_worker, n_rows)
    parity2t = generate_empty_matrix(n_cols_per_worker, n_rows)
    for i in range(1, n_procs-2):
        sub_matrix = feature_matrix[(i-1)*n_rows_per_worker:(i)*n_rows_per_worker]
        save_matrix(sub_matrix, output_dir + str(i) + ".dat")
        sub_matrix_transpose = transpose_matrix[(i-1)*n_cols_per_worker:(i)*n_cols_per_worker]
        save_matrix(sub_matrix_transpose, output_dir + str(i) + "t.dat")

        sum_matrices(parity1, sub_matrix, n_rows_per_worker, n_cols)
        multmat = scal_mult_matrix_ret(sub_matrix, i, n_rows_per_worker, n_cols)
        sum_matrices(parity2, multmat, n_rows_per_worker, n_cols)

        sum_matrices(parity1t, sub_matrix_transpose, n_cols_per_worker, n_rows)
        multmatt = scal_mult_matrix_ret(sub_matrix_transpose, i, n_cols_per_worker, n_rows)
        sum_matrices(parity2t, multmatt, n_cols_per_worker, n_rows)
    save_matrix(parity1, output_dir + str(n_procs-2) + ".dat")
    save_matrix(parity2, output_dir + str(n_procs-1) + ".dat")
    save_matrix(parity1t, output_dir + str(n_procs-2) + "t.dat")
    save_matrix(parity2t, output_dir + str(n_procs-1) + "t.dat")
    save_vector(labels_vector, output_dir + "labels.dat")
    regr = linear_model.LinearRegression()
    regr.fit(feature_matrix, labels_vector)
    answer = regr.coef_
    save_vector(answer, output_dir + "answer.dat")
