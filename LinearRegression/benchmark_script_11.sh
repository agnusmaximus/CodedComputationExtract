N_PROCS=11
ROWS=(2000 400 10000)
COLS=(2000 10000 400)
OUTPUT_PATH=../../data/linear_regression/m1small/
OUTPUT_DIR=${OUTPUT_PATH}${N_PROCS}_processes/
N_ITERS=100

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

for i in $(seq 0 $((${#COLS[*]}-1)))
do
    echo "Running (row) 2-coded linear regression on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded2_${ROWS[i]}_${COLS[i]}
    make N_PROCS=${N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} M_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} run_coded_many

    echo "Running uncoded linearregression on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}rowuncoded_${ROWS[i]}_${COLS[i]}
    make N_PROCS=${N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} M_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} run_uncoded_many
done
