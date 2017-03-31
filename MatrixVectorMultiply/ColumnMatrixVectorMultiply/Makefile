N_ROWS=720
N_COLS=720
N_PROCS=11
OUTPUT_FILE=uncoded_column
N_COMPARE_RUNS=10

uncoded_run:
	python data/generate_col_partitioned_matrix.py $(N_PROCS) 0 $(N_COLS) $(N_ROWS) data/$(N_ROWS)x$(N_COLS)uncoded_column 1
	mpic++ ColumnMatrixVectorMultiply.cpp -DCODED=0 -DDATA_PATH=\"data/$(N_ROWS)x$(N_COLS)uncoded_column/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lcblas -o uncoded_column
	rm -rf ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./uncoded_column >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}
