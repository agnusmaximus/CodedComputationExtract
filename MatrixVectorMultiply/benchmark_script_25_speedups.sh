# ROW UNCODED PARAMS
N_PROCS_ROW=(5 7 9 13 17 21 23 25)
ROWS=(5748 5748 5744 5748 5744 5740 5742 5736)
COLS=(5748 5748 5744 5748 5744 5740 5742 5736)

# CODED 1 PARAMS
N_PROCS_CODED_1=(5 7 9 13 17 21 23 25)
ROWS_CODED_1=(5751 5750 5754 5753 5760 5757 5754 5750)
COLS_CODED_1=(5751 5750 5754 5753 5760 5757 5754 5750)

# CODED 2 PARAMS
N_PROCS_CODED_2=(5 7 9 13 17 21 23 25)
ROWS_CODED_2=(5750 5752 5754 5750 5754 5760 5760 5764)
COLS_CODED_2=(5750 5752 5754 5750 5754 5760 5760 5764)

#OUTPUT_PATH=/home/ubuntu/CodedComputation/data/c1medium/
OUTPUT_PATH=/Users/maxlam/Desktop/
OUTPUT_DIR=${OUTPUT_PATH}${N_PROCS}_processes/
N_ITERS=100

rm -rf ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}

for i in $(seq 0 $((${#N_PROCS_ROW[*]}-1)))
do

     echo "Running row uncoded matrix vector multiplication on " ${ROWS[i]} ${COLS[i]} " matricse with " ${N_PROCS_ROW[i]} " workers "
     OUTPUT_FILE=${OUTPUT_DIR}row_uncoded_${ROWS[i]}_${COLS[i]}_nprocs_${N_PROCS_ROW[i]}
     cd RowMatrixVectorMultiply/
     make N_PROCS=${N_PROCS_ROW[i]} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS[i]} N_COLS=${COLS[i]} OUTPUT_FILE=${OUTPUT_FILE} uncoded_run
     cd data
     find -not -iname "*.py" | xargs rm -rf
     cd ../
     cd ../


    echo "Running coded-1 matrix vector multiplication on " ${ROWS_CODED_1[i]} ${COLS_CODED_1[i]} " matricse with " ${N_PROCS_CODED_1[i]} " workers "
    OUTPUT_FILE=${OUTPUT_DIR}coded_1_${ROWS_CODED_1[i]}_${COLS_CODED_1[i]}_nprocs_${N_PROCS_CODED_1[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_CODED_1[i]} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_CODED_1[i]} N_COLS=${COLS_CODED_1[i]} OUTPUT_FILE=${OUTPUT_FILE} coded1_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

    echo "Running coded-2 matrix vector multiplication on " ${ROWS_CODED_2[i]} ${COLS_CODED_2[i]} " matricse with " ${N_PROCS_CODED_2[i]} " workers "
    OUTPUT_FILE=${OUTPUT_DIR}coded_2_${ROWS_CODED_2[i]}_${COLS_CODED_2[i]}_nprocs_${N_PROCS_CODED_2[i]}
    cd RowMatrixVectorMultiply/
    make N_PROCS=${N_PROCS_CODED_2[i]} N_COMPARE_RUNS=${N_ITERS} N_ROWS=${ROWS_CODED_2[i]} N_COLS=${COLS_CODED_2[i]} OUTPUT_FILE=${OUTPUT_FILE} coded2_run
    cd data
    find -not -iname "*.py" | xargs rm -rf
    cd ../
    cd ../

done
