import sys
import numpy as np
import scipy.io as io

def load_input_file(fname):
    f = open(fname)
    arr = []
    for line in f:
        arr.append(int(line))
    print(arr)
    f.close()
    return arr

input_file, output_file_name = sys.argv[1], sys.argv[2]
arr = np.array(load_input_file(input_file))
d = {"data":arr}
io.savemat(output_file_name + ".mat",d)
