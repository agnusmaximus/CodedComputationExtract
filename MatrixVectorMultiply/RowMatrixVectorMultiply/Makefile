#N_ROWS=2300
#N_COLS=2300
N_ROWS=960
N_COLS=960
N_PROCS=25
N_COMPARE_RUNS=50
NUM_RUNS=1
OUTPUT_PATH=../../data/worker_runs/$(N_ROWS)x$(N_COLS)/
OUTPUT_FILE=row_mult_output

# For replication
N_CHUNKS = 5
N_WORKERS_PER_CHUNK = 5

# For simulated parity
PARITY = 3

# For lrc
N_SPLITS = 6

lrc_coded_run:
	python data/generate_partitioned_matrix_lrc.py $(N_PROCS) $(N_SPLITS) $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)lrc$(N_SPLITS) 1 1
	mpic++ MatrixVectorMultiply.cpp -DN_SPLITS=$(N_SPLITS) -DCODED=-5 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)lrc$(N_SPLITS)/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o lrc_coded
	rm -rf ${OUTPUT_FILE}
	#mpirun -np $(N_PROCS) ./lrc
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./lrc_coded >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}

lrc_run:
	python data/generate_partitioned_matrix_lrc.py $(N_PROCS) $(N_SPLITS) $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)lrc$(N_SPLITS) 1 0
	mpic++ MatrixVectorMultiply.cpp -DN_SPLITS=$(N_SPLITS) -DCODED=-4 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)lrc$(N_SPLITS)/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o lrc
	rm -rf ${OUTPUT_FILE}
	#mpirun -np $(N_PROCS) ./lrc
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./lrc >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}

replicate_run:
	python data/generate_partitioned_matrix_replicate.py $(N_PROCS) $(N_CHUNKS) $(N_WORKERS_PER_CHUNK) $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)replicate$(N_CHUNKS)x$(N_WORKERS_PER_CHUNK) 1
	mpic++ MatrixVectorMultiply.cpp -DN_CHUNKS=$(N_CHUNKS) -DN_WORKERS_PER_CHUNK=$(N_WORKERS_PER_CHUNK) -DCODED=-2 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)replicate$(N_CHUNKS)x$(N_WORKERS_PER_CHUNK)/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o replicated
	rm -rf ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./replicated >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}
simulated_coded_run:
	python data/generate_partitioned_matrix_simulate.py $(N_PROCS) $(PARITY) $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)simulated_coded$(PARITY) 0
	mpic++ MatrixVectorMultiply.cpp -DCODED=-3 -DPARITY=$(PARITY) -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)simulated_coded$(PARITY)/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o simulated_coded
	rm -rf ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./simulated_coded >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}
coded2_run:
	python data/generate_partitioned_matrix.py $(N_PROCS) 2 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded2 1
	mpic++ MatrixVectorMultiply.cpp -DCODED=2 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded2/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded2
	rm -rf ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded2 >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}

coded1_run:
	python data/generate_partitioned_matrix.py $(N_PROCS) 1 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded 1
	mpic++ MatrixVectorMultiply.cpp -DCODED=1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded1
	rm -rf ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded1 >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}

uncoded_run:
	python data/generate_partitioned_matrix.py $(N_PROCS) 0 $(N_COLS) $(N_ROWS) ./data/$(N_ROWS)x$(N_COLS)uncoded 1
	mpic++ MatrixVectorMultiply.cpp -DNUM_RUNS=$(NUM_RUNS) -DCODED=0 -DDATA_PATH=\"./data/$(N_ROWS)x$(N_COLS)uncoded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lblas -o uncoded
	rm -rf ${OUTPUT_FILE}
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./uncoded >> ${OUTPUT_FILE}; done
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'


	rm -rf tmp_out
	more ${OUTPUT_FILE} | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' > tmp_out
	python ../save_as_mat.py tmp_out ${OUTPUT_FILE}

all:
	mpic++ MatrixVectorMultiply.cpp -DNUM_RUNS=$(NUM_RUNS) -DCODED=0 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)uncoded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lblas -o uncoded
	mpic++ MatrixVectorMultiply.cpp -DNUM_RUNS=$(NUM_RUNS) -DCODED=1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lblas -o coded
	mpic++ MatrixVectorMultiply.cpp -DNUM_RUNS=$(NUM_RUNS) -DCODED=2 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded2/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lblas -o coded2
	#mpic++ MatrixVectorMultiply.cpp -DCODED=-1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)extreme/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lblas -o extreme
compare:
	python data/generate_partitioned_matrix.py $(N_PROCS) 0 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)uncoded 1
	python data/generate_partitioned_matrix.py $(N_PROCS) 1 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded 1
	python data/generate_partitioned_matrix.py $(N_PROCS) 2 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded2 1
	#python data/generate_partitioned_matrix.py $(N_PROCS) -1 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)extreme 1
	mpic++ MatrixVectorMultiply.cpp -DCODED=0 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)uncoded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o uncoded
	mpic++ MatrixVectorMultiply.cpp -DCODED=1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded
	mpic++ MatrixVectorMultiply.cpp -DCODED=2 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded2/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded2
	#mpic++ MatrixVectorMultiply.cpp -DCODED=-1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)extreme/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -I /usr/local/opt/openblas/include/ -lblas -o extreme
	mkdir -p $(OUTPUT_PATH)
	rm -f $(OUTPUT_PATH)uncoded_output
	rm -f $(OUTPUT_PATH)coded_output
	rm -f $(OUTPUT_PATH)coded2_output
	#rm -f extreme_output
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded >> $(OUTPUT_PATH)coded_output && mpirun -np $(N_PROCS) ./uncoded >> $(OUTPUT_PATH)uncoded_output && mpirun -np $(N_PROCS) ./coded2 >> $(OUTPUT_PATH)coded2_output; done
	#for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./uncoded >> $(OUTPUT_PATH)uncoded_output; done
	#for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded2 >> $(OUTPUT_PATH)coded2_output; done
	#for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./extreme >> extreme_output; done
	more $(OUTPUT_PATH)coded2_output | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	more $(OUTPUT_PATH)coded_output | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	more $(OUTPUT_PATH)uncoded_output | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	#more extreme_output | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	more coded2_output | grep "WRONG" | head
	more coded_output | grep "WRONG" | head
	more uncoded_output | grep "WRONG" | head
	#more extreme_output | grep "WRONG" | head
compare_runs:
	python data/generate_partitioned_matrix.py $(N_PROCS) 0 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)uncoded 1
	python data/generate_partitioned_matrix.py $(N_PROCS) 2 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded2 1
	mpic++ MatrixVectorMultiply.cpp -DCODED=3 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded2/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded2_runs
	mpic++ MatrixVectorMultiply.cpp -DCODED=4 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)uncoded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS)  -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o uncoded_runs
	rm -f uncoded_output_runs$(NUM_RUNS)
	rm -f coded2_output_runs$(NUM_RUNS)
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./uncoded_runs >> uncoded_output_runs$(NUM_RUNS); done
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded2_runs >> coded2_output_runs$(NUM_RUNS); done
	more coded2_output_runs$(NUM_RUNS) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	more uncoded_output_runs$(NUM_RUNS) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
profile_coded:
	python data/generate_partitioned_matrix.py $(N_PROCS) 1 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded 1
	mpic++ MatrixVectorMultiply.cpp -DPROFILE=1 -DCODED=1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded_profile
	mpirun -np $(N_PROCS) ./coded_profile
profile_coded2:
	python data/generate_partitioned_matrix.py $(N_PROCS) 2 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded2 1
	mpic++ MatrixVectorMultiply.cpp -DPROFILE=1 -DCODED=2 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded2/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded2_profile
	mpirun -np $(N_PROCS) ./coded2_profile
profile_uncoded:
	python data/generate_partitioned_matrix.py $(N_PROCS) 0 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)uncoded 1
	mpic++ MatrixVectorMultiply.cpp -DPROFILE=1 -DCODED=0 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)uncoded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o uncoded_profile
	mpirun -np $(N_PROCS) ./uncoded_profile
profile_all:
	python data/generate_partitioned_matrix.py $(N_PROCS) 1 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded 1
	mpic++ MatrixVectorMultiply.cpp -DPROFILE=1 -DCODED=1 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded_profile
	python data/generate_partitioned_matrix.py $(N_PROCS) 2 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)coded2 1
	mpic++ MatrixVectorMultiply.cpp -DPROFILE=1 -DCODED=2 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)coded2/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o coded2_profile
	python data/generate_partitioned_matrix.py $(N_PROCS) 0 $(N_COLS) $(N_ROWS) /data/$(N_ROWS)x$(N_COLS)uncoded 1
	mpic++ MatrixVectorMultiply.cpp -DPROFILE=1 -DCODED=0 -DDATA_PATH=\"/data/$(N_ROWS)x$(N_COLS)uncoded/\" -DN_ROWS=$(N_ROWS) -DN_COLS=$(N_COLS) -DNUM_RUNS=$(NUM_RUNS) -I /usr/local/opt/openblas/include/ -lblas -o uncoded_profile

	rm -f profile_coded
	rm -f profile_coded2
	rm -f profile_uncoded

	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded_profile >> profile_coded; done
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./coded2_profile >> profile_coded2; done
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./uncoded_profile >> profile_uncoded; done

	@echo 'CODED 2'
	@echo '----------------------------'
	@echo 'RECEIVE INPUT TIME: '
	@more profile_coded2 | grep "RECEIVE INPUT TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	@echo 'WORK TIME: '
	@more profile_coded2 | grep "WORK TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	@echo 'SEND RESULT TIME: '
	@more profile_coded2 | grep "SEND RESULT TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	@echo 'CODED 1'
	@echo '----------------------------'
	@echo 'RECEIVE INPUT TIME: '
	@more profile_coded | grep "RECEIVE INPUT TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	@echo 'WORK TIME: '
	@more profile_coded | grep "WORK TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	@echo 'SEND RESULT TIME: '
	@more profile_coded | grep "SEND RESULT TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'

	@echo 'UNCODED'
	@echo '----------------------------'
	@echo 'RECEIVE INPUT TIME: '
	@more profile_uncoded | grep "RECEIVE INPUT TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	@echo 'WORK TIME: '
	@more profile_uncoded | grep "WORK TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	@echo 'SEND RESULT TIME: '
	@more profile_uncoded | grep "SEND RESULT TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
