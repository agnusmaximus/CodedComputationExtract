N_ROWS=5750
N_COLS=5750
N_COMPARE_RUNS=100

OUTPUT_FILE=row_mult_output

single_machine_run:
	python data/generate_partitioned_matrix.py $(N_ROWS) $(N_COLS) data/$(N_ROWS)x$(N_COLS)single_machine
	g++ SingleMachineMatrixVectorMultiply.cpp -DDATA_PATH=\"data/$(N_ROWS)x$(N_COLS)single_machine/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -I /usr/local/opt/openblas/include/ -lblas -o single_machine
	rm -f ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do  ./single_machine >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}


single_machine_run_single:
	python data/generate_partitioned_matrix.py $(N_ROWS) $(N_COLS) data/$(N_ROWS)x$(N_COLS)single_machine
	g++ SingleMachineMatrixVectorMultiply.cpp -DDATA_PATH=\"data/$(N_ROWS)x$(N_COLS)single_machine/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -I /usr/local/opt/openblas/include/ -lblas -o single_machine
	./single_machine
