# UNCODED PARAMS
N_PROCS=25

ROWS=(21240)
COLS=(21240)

# REAL CODED PARAMS
CODED_2_ROWS=(21274)
CODED_2_COLS=(21274)
CODED_2_N_PROCS=25

CODED_1_ROWS=(21275)
CODED_1_COLS=(21275)
CODED_1_N_PROCS=25

LRC_6_ROWS=(21264)
LRC_6_COLS=(21264)
LRC_6_N_PROCS=25

# FAKE CODED PARAMS
FAKE_CODED_2_ROWS=(21274)
FAKE_CODED_2_COLS=(21274)
FAKE_CODED_2_N_PROCS=25
FAKE_CODED_2_PARITY=2

FAKE_CODED_4_ROWS=(21260)
FAKE_CODED_4_COLS=(21260)
FAKE_CODED_4_N_PROCS=25
FAKE_CODED_4_PARITY=4

FAKE_CODED_6_ROWS=(21276)
FAKE_CODED_6_COLS=(21276)
FAKE_CODED_6_N_PROCS=25
FAKE_CODED_6_PARITY=6

FAKE_CODED_8_ROWS=(21264)
FAKE_CODED_8_COLS=(21264)
FAKE_CODED_8_N_PROCS=25
FAKE_CODED_8_PARITY=8

FAKE_CODED_10_ROWS=(21266)
FAKE_CODED_10_COLS=(21266)
FAKE_CODED_10_N_PROCS=25
FAKE_CODED_10_PARITY=10

FAKE_CODED_12_ROWS=(5760)
FAKE_CODED_12_COLS=(5760)
FAKE_CODED_12_N_PROCS=25
FAKE_CODED_12_PARITY=12

FAKE_CODED_14_ROWS=(21260)
FAKE_CODED_14_COLS=(21260)
FAKE_CODED_14_N_PROCS=25
FAKE_CODED_14_PARITY=14

FAKE_CODED_16_ROWS=(21264)
FAKE_CODED_16_COLS=(21264)
FAKE_CODED_16_N_PROCS=25
FAKE_CODED_16_PARITY=16

# REPLICATION PARAMS
ROWS_REPL_12_2=(21252)
COLS_REPL_12_2=(21252)
N_PROCS_12_2=25

ROWS_REPL_8_3=(21256)
COLS_REPL_8_3=(21256)
N_PROCS_8_3=25

ROWS_REPL_6_4=(21258)
COLS_REPL_6_4=(21258)
N_PROCS_6_4=25

ROWS_REPL_4_6=(21260)
COLS_REPL_4_6=(21260)
N_PROCS_4_6=25

ROWS_REPL_3_8=(21258)
COLS_REPL_3_8=(21258)
N_PROCS_3_8=25

#OUTPUT_PATH=/home/ubuntu/CodedComputation/data/c1medium/
OUTPUT_PATH=/Users/maxlam/Desktop/
OUTPUT_DIR=${OUTPUT_PATH}${N_PROCS}_processes/
N_ITERS=50

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

for i in $(seq 0 $((${#COLS[*]}-1)))
do

    ############################
    # REAL CODED               #
    ############################
    echo "Running coded lrc 6-split (~12coded) matrix vector multiplication on " ${LRC_6_ROWS[i]} ${LRC_6_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded_lrc6_${LRC_6_ROWS[i]}_${LRC_6_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=6 N_PROCS=${LRC_6_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${LRC_6_ROWS[i]} N_COLS=${LRC_6_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} lrc_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running (row) 1-coded matrix vector multiplication on" ${CODED_1_ROWS[i]} ${CODED_1_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded1_${CODED_1_ROWS[i]}_${CODED_1_COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${CODED_1_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${CODED_1_ROWS[i]} N_COLS=${CODED_1_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} coded1_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running (row) 2-coded matrix vector multiplication on" ${CODED_2_ROWS[i]} ${CODED_2_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}coded2_${CODED_2_ROWS[i]}_${CODED_2_COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${CODED_2_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${CODED_2_ROWS[i]} N_COLS=${CODED_2_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} coded2_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    ##############################
    # FAKE CODED                 #
    ##############################
    echo "Running fake 2-coded matrix vector multiplication on" ${FAKE_CODED_2_ROWS[i]} ${FAKE_CODED_2_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded2_${FAKE_CODED_2_ROWS[i]}_${FAKE_CODED_2_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_2_PARITY} N_PROCS=${FAKE_CODED_2_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_2_ROWS[i]} N_COLS=${FAKE_CODED_2_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 4-coded matrix vector multiplication on" ${FAKE_CODED_4_ROWS[i]} ${FAKE_CODED_4_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded4_${FAKE_CODED_4_ROWS[i]}_${FAKE_CODED_4_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_4_PARITY} N_PROCS=${FAKE_CODED_4_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_4_ROWS[i]} N_COLS=${FAKE_CODED_4_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 6-coded matrix vector multiplication on" ${FAKE_CODED_6_ROWS[i]} ${FAKE_CODED_6_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded6_${FAKE_CODED_6_ROWS[i]}_${FAKE_CODED_6_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_6_PARITY} N_PROCS=${FAKE_CODED_6_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_6_ROWS[i]} N_COLS=${FAKE_CODED_6_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 8-coded matrix vector multiplication on" ${FAKE_CODED_8_ROWS[i]} ${FAKE_CODED_8_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded8_${FAKE_CODED_8_ROWS[i]}_${FAKE_CODED_8_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_8_PARITY} N_PROCS=${FAKE_CODED_8_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_8_ROWS[i]} N_COLS=${FAKE_CODED_8_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 10-coded matrix vector multiplication on" ${FAKE_CODED_10_ROWS[i]} ${FAKE_CODED_10_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded10_${FAKE_CODED_10_ROWS[i]}_${FAKE_CODED_10_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_10_PARITY} N_PROCS=${FAKE_CODED_10_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_10_ROWS[i]} N_COLS=${FAKE_CODED_10_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 12-coded matrix vector multiplication on" ${FAKE_CODED_12_ROWS[i]} ${FAKE_CODED_12_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded12_${FAKE_CODED_12_ROWS[i]}_${FAKE_CODED_12_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_12_PARITY} N_PROCS=${FAKE_CODED_12_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_12_ROWS[i]} N_COLS=${FAKE_CODED_12_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 14-coded matrix vector multiplication on" ${FAKE_CODED_14_ROWS[i]} ${FAKE_CODED_14_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded14_${FAKE_CODED_14_ROWS[i]}_${FAKE_CODED_14_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_14_PARITY} N_PROCS=${FAKE_CODED_14_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_14_ROWS[i]} N_COLS=${FAKE_CODED_14_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running fake 16-coded matrix vector multiplication on" ${FAKE_CODED_16_ROWS[i]} ${FAKE_CODED_16_COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}fake_coded16_${FAKE_CODED_16_ROWS[i]}_${FAKE_CODED_16_COLS[i]}
    cd RowMatrixVectorMultiply/
    make PARITY=${FAKE_CODED_16_PARITY} N_PROCS=${FAKE_CODED_16_N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${FAKE_CODED_16_ROWS[i]} N_COLS=${FAKE_CODED_16_COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} simulated_coded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    ###############################
    # UNCODED                     #
    ###############################
    echo "Running row matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}rowuncoded_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running column matrix vector multiplication on" ${ROWS[i]} ${COLS[i]} "matrices"
    OUTPUT_FILE=${OUTPUT_DIR}columnuncoded_${ROWS[i]}_${COLS[i]}
    cd ColumnMatrixVectorMultiply/
    make N_PROCS=${N_PROCS} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
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

    ###############################
    # REPLICATED                  #
    ##############################
    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_12_2[i]} ${COLS_REPL_12_2[i]} " matrices with 12_2 (chunks, workers) per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_12_2_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=12 N_WORKERS_PER_CHUNK=2 N_PROCS=${N_PROCS_12_2} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_12_2[i]} N_COLS=${COLS_REPL_12_2[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_8_3[i]} ${COLS_REPL_8_3[i]} " matrices with 8_3 (chunks, workers) per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_8_3_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=8 N_WORKERS_PER_CHUNK=3 N_PROCS=${N_PROCS_8_3} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_8_3[i]} N_COLS=${COLS_REPL_8_3[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_6_4[i]} ${COLS_REPL_6_4[i]} " matrices with 6_4 (chunks, workers) per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_6_4_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=6 N_WORKERS_PER_CHUNK=4 N_PROCS=${N_PROCS_6_4} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_6_4[i]} N_COLS=${COLS_REPL_6_4[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_4_6[i]} ${COLS_REPL_4_6[i]} " matrices with 4_6 (chunks, workers) per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_4_6_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=4 N_WORKERS_PER_CHUNK=6 N_PROCS=${N_PROCS_4_6} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_4_6[i]} N_COLS=${COLS_REPL_4_6[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running row replicated uncoded matrix vector multiplication on" ${ROWS_REPL_3_8[i]} ${COLS_REPL_3_8[i]} " matrices with 3_8 (chunks, workers) per chunk"
    OUTPUT_FILE=${OUTPUT_DIR}rep_3_8_${ROWS[i]}_${COLS[i]}
    cd RowMatrixVectorMultiply/
    make N_CHUNKS=3 N_WORKERS_PER_CHUNK=8 N_PROCS=${N_PROCS_3_8} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_REPL_3_8[i]} N_COLS=${COLS_REPL_3_8[i]} OUTPUT_FILE=${OUTPUT_FILE} replicate_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../


done
