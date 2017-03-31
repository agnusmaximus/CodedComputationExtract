N_PROCS=36
N_COLS=36
N_ROWS=36
N_COMPARE_RUNS=10
OUTPUT_FILE = ./blocked_output

uncoded_run:
	python data/generate_partitioned_matrix.py $(N_ROWS) $(N_COLS)
	mpic++ -DDATA_PATH=\"\" mv_mult_mult_checkerboard.cpp  -I /usr/local/opt/openblas/include/ -lblas -o blocked
	rm -rf $(OUTPUT_FILE)
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./blocked >> $(OUTPUT_FILE); done
	more $(OUTPUT_FILE) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}
