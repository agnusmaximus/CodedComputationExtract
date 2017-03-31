from __future__ import print_function
import sys
import random
import os
import copy

def column_append(a, b):
    return [x + b[i] for i,x in enumerate(a)]

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

def generate_random_vector(n):
    return [random.randint(0, 100) for x in range(n)]


def generate_random_matrix(n_rows, n_cols):
    return [[random.randint(0, 100) for i in range(n_cols)] for x in range(n_rows)]

def matrix_vector_multiply(matrix, vector, n_rows, n_cols):
    result = [0 for x in range(n_rows)]
    for i in range(n_rows):
        for j in range(n_cols):
            result[i] += matrix[i][j] * vector[j]
    return result

if len(sys.argv) != 7:
    print("Usage: python generate_partitioned_matrix.py n_procs is_coded n_cols n_rows output_dir should_output_result")
    sys.exit(0)

n_procs, is_coded, n_cols, n_rows, output_dir, should_calculate_answer  = [x for x in sys.argv[1:]]
n_procs, is_coded, n_cols, n_rows, should_calculate_answer = int(n_procs), int(is_coded), int(n_cols), int(n_rows), int(should_calculate_answer)
output_dir = output_dir + "/" if output_dir[-1] != "/" else output_dir

# num workers, and number of rows per worker
n_workers = n_procs - 1
n_cols_per_worker = n_cols / n_workers

# if generating coded partitioned matrix, make sure the number of rows is a multiple
# of n_workers-1
if is_coded:
    if is_coded == 1:
        assert(n_rows % (n_workers-1) == 0)
    else:
        assert(n_rows % (n_workers-2) == 0)
else:
    assert(n_rows % (n_workers)== 0)

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Total matrix, if should calculate answer
total = []

print("GENERATING PARTITIONED MATRIX OF SIZE %d x %d FOR %d workers" % (n_rows, n_cols, n_workers))
print("- EACH WORKER GETS %d x %d doubles, %f MBYTES EACH" % (n_rows, n_cols_per_worker, (n_cols_per_worker * n_rows * 8) / 1000000.0))

if not is_coded:
    for i in range(1, n_procs):
        random_matrix = generate_random_matrix(n_rows, n_cols_per_worker)
        if total == []:
            total = copy.deepcopy(random_matrix)
        else:
            total = column_append(total, random_matrix)
        save_matrix(random_matrix, output_dir + str(i) + ".dat")
random_vector = generate_random_vector(n_cols)
save_vector(random_vector, output_dir + "vec.dat")
print("DONE")
if should_calculate_answer:
    result = matrix_vector_multiply(total, random_vector, n_rows, n_cols)
    save_vector(result, output_dir + "result.dat")
