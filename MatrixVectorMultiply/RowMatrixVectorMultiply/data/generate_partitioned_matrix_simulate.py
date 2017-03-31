from __future__ import print_function
import sys
import random
import os
from generate_partitioned_matrix_helpers import *

if len(sys.argv) != 7:
    print("Usage: python generate_partitioned_matrix_replicate.py n_procs n_replications n_workers_per_replicate n_cols n_rows output_dir should_calculate_answer")
    sys.exit(0)

n_procs, coded_parity_val, n_cols, n_rows, output_dir, should_calculate_answer  = [x for x in sys.argv[1:]]
n_procs, coded_parity_val, n_cols, n_rows, should_calculate_answer = int(n_procs), int(coded_parity_val), int(n_cols), int(n_rows), int(should_calculate_answer)
output_dir = output_dir + "/" if output_dir[-1] != "/" else output_dir

# num workers, and number of rows per worker
n_workers = n_procs - 1
n_rows_per_worker = n_rows / (n_workers-coded_parity_val)

# if generating coded partitioned matrix, make sure the number of rows is a multiple
# of n_workers-1
assert(n_rows % (n_workers-coded_parity_val) == 0)

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Total matrix, if should calculate answer
total = []

print("GENERATING PARTITIONED MATRIX OF SIZE %d x %d FOR %d workers" % (n_rows, n_cols, n_workers))
print("- EACH WORKER GETS %d x %d doubles, %f MBYTES EACH" % (n_rows_per_worker, n_cols, (n_rows_per_worker * n_cols * 8) / 1000000.0))
random_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)

#for i in range(1, n_procs-coded_parity_val):
for i in range(1, n_procs):
    generate_random_matrix_inplace(random_matrix, n_rows_per_worker, n_cols)
    #total += random_matrix
    save_matrix(random_matrix, output_dir + str(i) + ".dat")

random_vector = generate_random_vector(n_cols)
save_vector(random_vector, output_dir + "vec.dat")
print("DONE")
