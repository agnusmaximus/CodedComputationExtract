from __future__ import print_function
import sys
import random
import os
from generate_partitioned_matrix_helpers import *


if len(sys.argv) != 8:
    print("Usage: python generate_partitioned_matrix_lrc.py n_procs n_splits n_cols n_rows output_dir should_output_result coded_plus (assume 2 coded per split) (coded_plus stands for whether the n_splits are coded too)")
    sys.exit(0)

# Assume 2 coded for every split
n_procs, n_splits, n_cols, n_rows, output_dir, should_calculate_answer, coded_plus  = [x for x in sys.argv[1:]]
n_procs, n_splits, n_cols, n_rows, should_calculate_answer, coded_plus = int(n_procs), int(n_splits), int(n_cols), int(n_rows), int(should_calculate_answer), int(coded_plus)
output_dir = output_dir + "/" if output_dir[-1] != "/" else output_dir

n_workers = n_procs - 1

n_splits_use = n_splits
if coded_plus:
    n_splits_use -= 2

# Make sure n_rows divisible by n_splits_use
assert(n_rows % n_splits_use == 0)
# Make sure n_workers is divisible by n_splits_use
assert(n_workers % n_splits == 0)

n_workers_per_split = n_workers / n_splits
sub_matrix_n_rows = n_rows / n_splits_use

# Make sure that for every sub matrix dimension, it is 2 codable (divisible by n_workers - 2)
assert(sub_matrix_n_rows % (n_workers_per_split-2) == 0)

n_rows_per_worker = sub_matrix_n_rows / (n_workers_per_split-2)

assert(n_rows_per_worker * (n_workers_per_split-2) * n_splits_use == n_rows)

# Total matrix if should calculate answer
print("GENERATING %d PARTITIONED MATRIX OF SIZE %d x %d FOR %d workers" % (n_splits, n_rows, n_cols, n_workers))
print("- EACH WORKER GETS %d x %d doubles, %f MBYTES EACH" % (n_rows_per_worker, n_cols, (n_rows_per_worker * n_cols * 8) / 1000000.0))

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

total = []
if not coded_plus:
    for split in range(n_splits):
        parity1_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
        parity2_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
        for proc in range(1, n_workers_per_split+1-2):
            proc_id = proc + split * n_workers_per_split
            random_matrix = generate_random_matrix(n_rows_per_worker, n_cols)
            total += random_matrix
            save_matrix(random_matrix, output_dir + str(proc_id) + ".dat")
            sum_matrices(parity1_matrix, random_matrix, n_rows_per_worker, n_cols)
            multmat = scal_mult_matrix_ret(random_matrix, proc, n_rows_per_worker, n_cols)
            sum_matrices(parity2_matrix, multmat, n_rows_per_worker, n_cols)
        save_matrix(parity1_matrix, output_dir + str((split+1)*n_workers_per_split-1) + ".dat")
        save_matrix(parity2_matrix, output_dir + str((split+1)*n_workers_per_split) + ".dat")
else:
    parity1_split = generate_empty_matrix(sub_matrix_n_rows, n_cols)
    parity2_split = generate_empty_matrix(sub_matrix_n_rows, n_cols)
    random_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
    # last 2 splits are coded
    for split in range(n_splits-2):
        split_matrix_total = []
        parity1_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
        parity2_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
        for proc in range(1, n_workers_per_split+1-2):
            proc_id = proc + split * n_workers_per_split
            random_matrix = generate_random_matrix(n_rows_per_worker, n_cols)
            #total += random_matrix
            save_matrix(random_matrix, output_dir + str(proc_id) + ".dat")
            sum_matrices(parity1_matrix, random_matrix, n_rows_per_worker, n_cols)
            multmat = scal_mult_matrix_ret(random_matrix, proc, n_rows_per_worker, n_cols)
            sum_matrices(parity2_matrix, multmat, n_rows_per_worker, n_cols)
            split_matrix_total += random_matrix
        save_matrix(parity1_matrix, output_dir + str((split+1)*n_workers_per_split-1) + ".dat")
        save_matrix(parity2_matrix, output_dir + str((split+1)*n_workers_per_split) + ".dat")
        sum_matrices(parity1_split, split_matrix_total, sub_matrix_n_rows, n_cols)
        multmat_split = scal_mult_matrix_ret(split_matrix_total, split+1, sub_matrix_n_rows, n_cols)
        sum_matrices(parity2_split, multmat_split, sub_matrix_n_rows, n_cols)

    # Save last 2 parity splits
    parity1_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
    parity2_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
    for proc in range(1, n_workers_per_split+1-2):
        proc_id = proc + (n_splits-2) * n_workers_per_split
        partial = parity1_split[(proc-1)*n_rows_per_worker:(proc)*n_rows_per_worker]
        save_matrix(partial, output_dir + str(proc_id) + ".dat")
        sum_matrices(parity1_matrix, partial, n_rows_per_worker, n_cols)
        scal_mult_matrix_sum(parity2_matrix, multmat, proc, n_rows_per_worker, n_cols)
    save_matrix(parity1_matrix, output_dir + str((n_splits-2+1)*n_workers_per_split-1) + ".dat")
    save_matrix(parity2_matrix, output_dir + str((n_splits-2+1)*n_workers_per_split) + ".dat")
    parity1_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
    parity2_matrix = generate_empty_matrix(n_rows_per_worker, n_cols)
    for proc in range(1, n_workers_per_split+1-2):
        proc_id = proc + (n_splits-1) * n_workers_per_split
        partial = parity2_split[(proc-1)*n_rows_per_worker:(proc)*n_rows_per_worker]
        save_matrix(partial, output_dir + str(proc_id) + ".dat")
        sum_matrices(parity1_matrix, partial, n_rows_per_worker, n_cols)
        scal_mult_matrix_sum(parity2_matrix, multmat, proc, n_rows_per_worker, n_cols)
    save_matrix(parity1_matrix, output_dir + str((n_splits-1+1)*n_workers_per_split-1) + ".dat")
    save_matrix(parity2_matrix, output_dir + str((n_splits-1+1)*n_workers_per_split) + ".dat")

random_vector = generate_random_vector(n_cols)
save_vector(random_vector, output_dir + "vec.dat")
print("DONE")
assert(len(total) == n_rows)
if should_calculate_answer and 0:
    result = matrix_vector_multiply(total, random_vector, n_rows, n_cols)
    save_vector(result, output_dir + "result.dat")
