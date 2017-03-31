from __future__ import print_function
import sys
import random
import os
from generate_partitioned_matrix_helpers import *

n_rows, n_cols, path = sys.argv[1:]
n_rows, n_cols = [int(x) for x in [n_rows, n_cols]]
path = path + "/" if path[-1] != "/" else path

mat = generate_random_matrix(n_rows, n_cols)
vec = generate_random_vector(n_cols)

# Create output directory
if not os.path.exists(path):
    os.makedirs(path)

save_matrix(mat, path + str(0) + ".dat")
save_vector(vec, path + "vec.dat")
