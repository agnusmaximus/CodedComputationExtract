from __future__ import print_function
import sys
import random
import os
from generate_partitioned_matrix_helpers import *

if len(sys.argv) != 8:
    print("Usage: python generate_partitioned_matrix_replicate.py n_procs n_replications n_workers_per_replicate n_cols n_rows output_dir should_calculate_answer")
    sys.exit(0)

n_procs, n_chunks, n_workers_per_chunk, n_cols, n_rows, output_dir, should_calculate_answer = [x for x in sys.argv[1:]]
n_procs, n_chunks, n_workers_per_chunk, n_cols, n_rows, should_calculate_answer = (int(x) for x in [n_procs, n_chunks, n_workers_per_chunk, n_cols, n_rows, should_calculate_answer])
output_dir = output_dir + "/" if output_dir[-1] != "/" else output_dir

n_workers = n_procs-1
# The number of workers per piece * number of pieces = # of total workers
assert(n_chunks * n_workers_per_chunk == n_workers)
# Number of rows must be evenly divided by the number of chunks it is split into
assert(n_rows % n_chunks == 0)

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Calculate number of rows per worker
n_rows_per_worker = n_rows / n_chunks

# Total matrix, if should calculate answer
total = []

print("GENERATING REPLICATED PARTITIONED MATRIX OF SIZE %d x %d FOR %d WORKERS" % (n_rows, n_cols, n_workers))
print("- THERE ARE %d CHUNKS; EACH CHUNK HANDLED BY %d WORKERS" % (n_chunks, n_workers_per_chunk))
print("- EACH WORKER GETS %d x %d DOUBLES, %f MBYTES EACH" % (n_rows_per_worker, n_cols, (n_rows_per_worker * n_cols * 8) / 1000000.0))

# Generate n_chunks matrices
worker_id = 1
random_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
for i in range(n_chunks):
    generate_random_matrix_inplace(random_matrix, n_rows_per_worker, n_cols)
    total += random_matrix
    for j in range(n_workers_per_chunk):
        save_matrix(random_matrix, output_dir + str(worker_id) + ".dat")
        worker_id += 1

random_vector = generate_random_vector(n_cols)
save_vector(random_vector, output_dir + "vec.dat")
print("DONE")

if should_calculate_answer and 0:
    result = matrix_vector_multiply(total, random_vector, n_rows, n_cols)
    save_vector(result, output_dir + "result.dat")
