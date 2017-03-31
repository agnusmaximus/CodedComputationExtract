
/****************************************************************
                 C-DAC Tech Workshop : hyPack-2013
                        Oct 15-18,2013

   Example 5.7         : mv_mult_mult_checkerboard.c

   Objective           : Matrix_Vector Multiplication
                         (Using Block CheckerBoard Partitioning)

   Input               : Process 0 reads files (mdata.inp) for Matrix
                         and (vdata.inp) for Vector

   Output              : Process 0 prints the result of Matrix_Vector
                         Multiplication

  Created     	    : MAY-2012


   Necessary Condition : Number of processors should be perfect
                         square and less than or equal to 8

***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <time.h>
#include <sys/time.h>
#include "../util.cpp"

int main(int argc, char** argv)
{

  /* .......Variables Initialisation ......*/
  MPI_Status status;
  MPI_Comm row_comm, col_comm;
  MPI_Group MPI_GROUP_WORLD;

  int Numprocs, MyRank, root_p, Root = 0;
  int irow, icol, iproc, jproc, index, Proc_Id;
  int NoofRows, NoofCols, NoofRows_Bloc, NoofCols_Bloc;
  int Bloc_MatrixSize, Bloc_VectorSize, VectorSize;
  int Local_Index, Global_Row_Index, Global_Col_Index;
  long long int t1, t2;

  double **Matrix, *Matrix_Array, *Bloc_Matrix, *Vector, *Bloc_Vector;
  double *FinalResult, *MyResult, *FinalVector;

  int *ranks, colsize, colrank, rowsize, rowrank, myrow;
  int		MatrixFileStatus = 1, VectorFileStatus = 1;

  FILE *fp;


  /* MPI Initialisation  */
  MPI_Init(&argc, &argv);
  MPI_Comm_group(MPI_COMM_WORLD,&MPI_GROUP_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
  root_p = sqrt((double)Numprocs);
  if(Numprocs != (root_p * root_p)) {
     if( MyRank == 0)
        printf("Error : Number of processors should be perfect square \n");
     MPI_Finalize();
     exit(-1);
  }

  if(MyRank == 0)
  {
    /* .......Read the Input file ......*/
    if ((fp = fopen ("./data/mdata.inp", "r")) == NULL){
		 MatrixFileStatus = 0;
    }

    if(MatrixFileStatus != 0) {
    	fscanf(fp, "%d %d\n", &NoofRows,&NoofCols);

 	NoofRows_Bloc = NoofRows/root_p;
 	NoofCols_Bloc = NoofCols/root_p;

    /* ...Allocate memory and read Matrix from file .......*/
    	Matrix = (double **)malloc(NoofRows*sizeof(double *));
 	for(irow=0 ;irow < NoofRows; irow++){
	 	Matrix[irow] = (double *)malloc(NoofCols*sizeof(double));
	 	for(icol=0; icol < NoofCols; icol++) {
	    	     fscanf(fp, "%lf", &Matrix[irow][icol]);
       	        }
	}
	fclose(fp);
	 }

    /* Read vector from input file */
    if ((fp = fopen ("./data/vdata.inp", "r")) == NULL){
		 VectorFileStatus = 0;
    }

    if(VectorFileStatus != 0) {
    	fscanf(fp, "%d\n", &VectorSize);

	Vector = (double*)malloc(VectorSize*sizeof(double));
    	for(index = 0; index < VectorSize; index++)
       	    fscanf(fp, "%lf", &Vector[index]);

    }

   }  /* end  of if myrank = 0 */

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Bcast (&MatrixFileStatus, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(MatrixFileStatus == 0) {
      if(MyRank == Root)
	 printf("Can't open input file for Matrix ..... \n");
      MPI_Finalize();
      exit(-1);
   }

   MPI_Bcast (&VectorFileStatus, 1, MPI_INT, Root, MPI_COMM_WORLD);
   if(VectorFileStatus == 0) {
      if(MyRank == Root)
	 printf("Can't open input file for Vector ..... \n");
      MPI_Finalize();
      exit(-1);
   }

   MPI_Bcast (&NoofRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (&NoofCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (&VectorSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

   if(NoofCols != VectorSize){
      if(MyRank == 0){
	 printf("Matrice and Vector Dimensions incompatible for Multiplication");
      }
      MPI_Finalize();
      exit(-1);
   }

   if(NoofRows % root_p != 0 || NoofCols % root_p != 0){
      if(MyRank == 0)
	 printf("Matrix can't be divided among processors equally");
      MPI_Finalize();
      exit(-1);
   }

   NoofRows_Bloc = NoofRows / root_p;
   NoofCols_Bloc = NoofCols / root_p;

   if(MyRank == 0){

      Matrix_Array = (double *)malloc(NoofRows*NoofCols*sizeof(double));
      /* Convert Matrix into 1-D array according to Checkerboard ......*/
      Local_Index = 0;
      for(iproc = 0; iproc < root_p; iproc++){
         for(jproc = 0; jproc < root_p; jproc++){
   	     Proc_Id = iproc * root_p + jproc;
	     for(irow = 0; irow < NoofRows_Bloc; irow++){
		 Global_Row_Index = iproc * NoofRows_Bloc + irow;
	         for (icol = 0; icol < NoofCols_Bloc; icol++){
		      Global_Col_Index = jproc * NoofCols_Bloc + icol;
	              Matrix_Array[Local_Index++] =
			       Matrix[Global_Row_Index][Global_Col_Index];
	         }
	     }
	 }
       }
   }


  /* Memory allocating for Bloc Matrix */
  Bloc_VectorSize = VectorSize / root_p;
  Bloc_MatrixSize = NoofRows_Bloc * NoofCols_Bloc;
  Bloc_Matrix = (double *) malloc (Bloc_MatrixSize * sizeof(double));
  MPI_Scatter(Matrix_Array, Bloc_MatrixSize, MPI_DOUBLE, Bloc_Matrix,
	      Bloc_MatrixSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  /* Creating groups of procesors row wise */
  myrow=MyRank/root_p;
  MPI_Comm_split(MPI_COMM_WORLD,myrow,MyRank,&row_comm);
  MPI_Comm_size(row_comm,&rowsize);
  MPI_Comm_rank(row_comm,&rowrank);

  /* Creating groups of procesors column wise */
  myrow=MyRank%root_p;
  MPI_Comm_split(MPI_COMM_WORLD,myrow,MyRank,&col_comm);
  MPI_Comm_size(col_comm,&colsize);
  MPI_Comm_rank(col_comm,&colrank);

  if (MyRank == 0) t1 = get_time();

  /* Scatter part of vector to all row master processors */
  Bloc_Vector = (double*) malloc(Bloc_VectorSize * sizeof(double));
  if(MyRank/root_p == 0){
     MPI_Scatter(Vector, Bloc_VectorSize, MPI_DOUBLE, Bloc_Vector,
		 Bloc_VectorSize, MPI_DOUBLE, 0,row_comm);
  }

  /* Row master broadcasts its vector part to processors in its column */
  MPI_Bcast(Bloc_Vector, Bloc_VectorSize, MPI_DOUBLE, 0, col_comm);

  /* Multiplication done by all procs */

  MyResult   = (double*) malloc(NoofRows_Bloc * sizeof(double));
  index = 0;
  /*for(irow=0; irow < NoofRows_Bloc; irow++){
      MyResult[irow]=0;
      for(icol=0;icol< NoofCols_Bloc; icol++){
          MyResult[irow] += Bloc_Matrix[index++] * Bloc_Vector[icol];
      }
      }*/
  matrix_vector_multiply(Bloc_Matrix, Bloc_Vector, MyResult, NoofCols_Bloc, NoofRows_Bloc);


   /* collect partial product from all procs on to master processor
		and add it to get final answer */
   if(MyRank == Root)
      FinalResult = (double *)malloc(NoofRows_Bloc*Numprocs*sizeof(double));

   MPI_Gather (MyResult, NoofRows_Bloc, MPI_DOUBLE, FinalResult,
	       NoofRows_Bloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

/*   if(MyRank == 0){
   for(irow=0; irow < NoofRows_Bloc*Numprocs; irow++)
       printf(" * FinalVector[%d] = %7.3f\n",irow,FinalResult[irow]);
   }
*/

   if(MyRank == 0){
      FinalVector = (double *) malloc(NoofRows * sizeof(double));
      index = 0;
      for(iproc=0; iproc < root_p; iproc++){
	  for(irow=0; irow < NoofRows_Bloc; irow++){
	      FinalVector[index]  = 0;
	      for(jproc=0; jproc < root_p; jproc++){
		  FinalVector[index] +=
		    FinalResult[iproc*root_p*NoofRows_Bloc +
				jproc*NoofRows_Bloc +irow];

	      }
              //printf(" FinalVector[%d] = %7.3f\n",index,FinalVector[index]);
	      index++;
	  }
      }
      t2 = get_time();
      printf("ELAPSED TIME: %lld\n", t2-t1);
   }


   /* Free the groups formed */
   MPI_Comm_free(&row_comm);
   MPI_Comm_free(&col_comm);


   MPI_Finalize();

}
