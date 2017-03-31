N_PROCS=25
N_PROCS_ROW_COL=$(($N_PROCS+1))

ROWS=(5750 5750 11500)
COLS=(5750 11500 5750)

ROWS_REPL_13_2=(5720 5720 11440)
COLS_REPL_13_2=(5720 11440 5720)
N_PROCS_13_2=27

ROWS_REPL_9_3=(5670 5670 11340)
COLS_REPL_9_3=(5670 11340 5670)
N_PROCS_9_3=28

ROWS_REPL_7_4=(5740 5740 11480)
COLS_REPL_7_4=(5740 11480 5740)
N_PROCS_7_4=29

ROWS_REPL_5_5=(5750 5750 11500)
COLS_REPL_5_5=(5750 11500 5750)
N_PROCS_5_5=26

ROWS_SIMULATED_3=(5940 5940 11880)
COLS_SIMULATED_3=(5940 11880 5940)

ROWS_LRC_5=(5760 5760 11520)
COLS_LRC_5=(5760 11520 5760)

#OUTPUT_PATH=/home/ubuntu/CodedComputation/data/c1medium/
OUTPUT_PATH=/Users/maxlam/Desktop/
OUTPUT_DIR=${OUTPUT_PATH}${N_PROCS}_processes/
N_ITERS=100

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

for i in $(seq 0 $((${#COLS[*]}-1)))
do
    #echo "Running simulated 3-coded matrix vector multiplication on " ${ROWS_SIMULATED_3[i]} ${COLS_SIMULATED_3[i]} "matrices"
    #OUTPUT_FILE=${OUTPUT_DIR}simulated_coded3_${ROWS_SIMULATED_3[i]}_${COLS_SIMULATED_3[i]}
    #cd RowMatrixVectorMultiply/
    #make PARITY=3 N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_SIMULATED_3[i]} N_COLS=${COLS_SIMULATED_3[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    #cd data
    #find -not -iname "*.py" | xargs rm -rf
    #cd ../
    #cd ../

    echo "Running lrc 5-split (~10coded) matrix vector multiplication on " ${ROWS_LRC_5[i]} ${COLS_LRC_5[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}lrc5_${ROWS_LRC_5[i]}_${COLS_LRC_5[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=3 N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_LRC_5[i]} N_COLS=${COLS_LRC_5[i]} OUTPUT_FILE=${OUTPUT_FILE} lrc_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_13_2[i]} ${COLS_REPL_13_2[i]} " matrices with 13 chunks, 2 workers per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_13x2_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=13 N_WORKERS_PER_CHUNK=2 N_PROCS=${N_PROCS_13_2} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_13_2[i]} N_COLS=${COLS_REPL_13_2[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_9_3[i]} ${COLS_REPL_9_3[i]} " matrices with 9 chunks, 3 workers per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_9x3_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=9 N_WORKERS_PER_CHUNK=3 N_PROCS=${N_PROCS_9_3} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_9_3[i]} N_COLS=${COLS_REPL_9_3[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_7_4[i]} ${COLS_REPL_7_4[i]} " matrices with 7 chunks, 4 workers per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_7x4_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=7 N_WORKERS_PER_CHUNK=4 N_PROCS=${N_PROCS_7_4} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_7_4[i]} N_COLS=${COLS_REPL_7_4[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_5_5[i]} ${COLS_REPL_5_5[i]} " matrices with 5 chunks, 5 workers per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_5x5_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=5 N_WORKERS_PER_CHUNK=5 N_PROCS=${N_PROCS_5_5} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_5_5[i]} N_COLS=${COLS_REPL_5_5[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd ../

    echo "Running (row) 2-coded matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded2_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} coded2_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    #echo "Running (row) 1-coded matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    #OUTPUT_FILE=${OUTPUT_DIR}coded1_${ROWS[i]}_${COLS[i]}
    #cd RowMatrixVectorMultiply/
    #make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} coded1_run
    #cd ../

    echo "Running row matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}rowuncoded_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running column matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}columnuncoded_${ROWS[i]}_${COLS[i]}
    cd ColumnMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_ROW_COL} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running block matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}blockuncoded_${ROWS[i]}_${COLS[i]}
    cd BlockMatrixVectorMultiply
    make N_PROCS=${N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

done
