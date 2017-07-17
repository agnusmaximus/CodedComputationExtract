N_PROCS=11
M_COLS=2000
N_ROWS=2000
#M_COLS=400
#N_ROWS=10000
#M_COLS=10000
#N_ROWS=400
DATA_PATH=data
LEARNING_RATE= 2.46009504e-8 # Tuned for 2000x2000, 400x10000, 10000x400 mats
N_COMPARE_RUNS=20
N_ITERS=200
OUTPUT_FILE=out_file

run_coded_gantt:
	python data/generate_linear_regression_data.py ${N_PROCS} 2 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}coded" 1
	rm -rf gantt
	mkdir -p gantt
	mpic++ -DGANTT=1 -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=1 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}coded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_coded
	mpirun -np ${N_PROCS} ./lr_coded
	#python create_gantt_chart.py
run_coded_many:
	#python data/generate_linear_regression_data.py ${N_PROCS} 2 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}coded" 1
	mpic++ -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=1 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}coded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_coded
	rm -f $(OUTPUT_FILE)
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np ${N_PROCS} ./lr_coded >> $(OUTPUT_FILE); done
	more $(OUTPUT_FILE) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
run_coded_many_cancel:
	#python data/generate_linear_regression_data.py ${N_PROCS} 2 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}coded" 1
	mpic++ -DCANCEL_LG=1 -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=1 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}coded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_coded
	rm -f $(OUTPUT_FILE)
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np ${N_PROCS} ./lr_coded >> $(OUTPUT_FILE); done
	more $(OUTPUT_FILE) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
run_uncoded_many:
	python data/generate_linear_regression_data.py ${N_PROCS} 0 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}uncoded" 1
	mpic++ -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=0 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}uncoded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_uncoded
	rm -f $(OUTPUT_FILE)
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np ${N_PROCS} ./lr_uncoded >> $(OUTPUT_FILE); done
	more $(OUTPUT_FILE) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
all:
	mpic++ -DM_COLS=${M_COLS} -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr
run_uncoded:
	python data/generate_linear_regression_data.py ${N_PROCS} 0 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}uncoded" 1
	mpic++ -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=0 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}uncoded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_uncoded
	mpirun -np ${N_PROCS} ./lr_uncoded
run_coded:
	#python data/generate_linear_regression_data.py ${N_PROCS} 2 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}coded" 1
	mpic++ -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=1 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}coded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_coded
	mpirun -np ${N_PROCS} ./lr_coded
run_coded_cancellation:
	#python data/generate_linear_regression_data.py ${N_PROCS} 2 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}coded" 1
	mpic++ -DCANCEL_LG=1 -DN_ITERS=${N_ITERS} -DM_COLS=${M_COLS} -DTYPE=1 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}coded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_coded
	mpirun -np ${N_PROCS} ./lr_coded
compare:
	python data/generate_linear_regression_data.py ${N_PROCS} 0 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}uncoded" 1
	python data/generate_linear_regression_data.py ${N_PROCS} 2 ${M_COLS} ${N_ROWS} "${DATA_PATH}/${N_ROWS}x${M_COLS}coded" 1
	mpic++ -DM_COLS=${M_COLS} -DTYPE=1 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}coded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_coded
	mpic++ -DM_COLS=${M_COLS} -DTYPE=0 -DN_ROWS=${N_ROWS} -DDATA_PATH=\"${DATA_PATH}/${N_ROWS}x${M_COLS}uncoded/\" -DLEARNING_RATE=${LEARNING_RATE} LinearRegression.cpp  -I /usr/local/opt/openblas/include/ -lcblas -o lr_uncoded
	rm -f uncoded_runs$(N_COMPARE_RUNS)
	rm -f coded2_runs$(N_COMPARE_RUNS)
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./lr_uncoded >> uncoded_runs$(N_COMPARE_RUNS); done
	for i in `seq 1 $(N_COMPARE_RUNS)`; do mpirun -np $(N_PROCS) ./lr_coded >> coded2_runs$(N_COMPARE_RUNS); done
	more coded2_runs$(N_COMPARE_RUNS) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
	more uncoded_runs$(N_COMPARE_RUNS) | grep "ELAPSED TIME: " | grep -oEi '[0-9]+' | awk 'NR == 1 { max=$$1; min=$$1; sum=0 } { if ($$1>max) max=$$1; if ($$1<min) min=$$1; sum+=$$1;} END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}'
