N_PROCS=16
N_PROCS_ROW_COL=$(($N_PROCS+1))
ROWS=(3360 3360 6720)
COLS=(3360 6720 3360)
OUTPUT_PATH=/home/ubuntu/CodedComputation/data/c1medium/
OUTPUT_DIR=${OUTPUT_PATH}${N_PROCS}_processes/
N_ITERS=100

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

for i in $(seq 0 $((${#COLS[*]}-1)))
do
    echo "Running (row) 2-coded matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded2_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} coded2_run
    cd ../

    echo "Running (row) 1-coded matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded1_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} coded1_run
    cd ../

    echo "Running row matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}rowuncoded_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd ../

    echo "Running column matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}columnuncoded_${ROWS[i]}_${COLS[i]}
    cd ColumnMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd ../

    echo "Running block matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}blockuncoded_${ROWS[i]}_${COLS[i]}
    cd BlockMatrixVectorMultiply

    make N_PROCS=${N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd ../
done
