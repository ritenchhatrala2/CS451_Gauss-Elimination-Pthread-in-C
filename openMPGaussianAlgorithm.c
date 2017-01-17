/*
 * openMP.c
 *
 *  Created on: Oct 1, 2016
 *      Author: riten
 */
#include <stdio.h>      /* Input/Output */
#include <stdlib.h>     /* General Utilities */
#include <omp.h>		/* OpenMP Lib*/


/*************************************************
 * Structure to globalize the data for Thread Use.
 *************************************************/
struct Thread_Data {
	double** A;
	int n;
	int numThreads;
} t_info;

/*Declareration of Pthread and openMP function*/
void  openMPAlgorith( double** C );

/******************************************************
 * Generate random values for provided order mattrix
 ******************************************************/
void generateMatrix(int N, double** A){
	int i, j;
	for (i=0; i<N; i++){
		for (j=0; j<N+1; j++){
			A[i][j] = (rand()/50000);
		}
	}
}

/***********************************************************
 * print Matrix if required.large matrix will not be printed.
 ***********************************************************/
void printMarix(char* msg,double** A){
	int i,j;
	printf("%s\n",msg);
		for(i=0;i<t_info.n;i++)
		{
			for(j=0;j<t_info.n+1;j++)
			{
				printf("%f\t",A[i][j]);
			}
			printf("\n");
		}
}

/************************************************************
 * Back Substitution implementation after matrix transform to
 * upper triangular form
 ************************************************************/
void backSubstitution(double** A,double* x){
	int i,j;
	double sum;
	x[t_info.n-1]=A[t_info.n-1][t_info.n]/A[t_info.n-1][t_info.n-1];

		for(i=t_info.n-2; i>=0; i--)
		{
			sum=0;
			for(j=i+1; j<=t_info.n-1; j++)
			{
				sum=sum+A[i][j]*x[j];
			}
			x[i]=(A[i][t_info.n]-sum)/A[i][i];
		}
}

/*******************************************************************
 * Main() function execution start.
 * Serial, Pthread , Pthread Version 2 and openMP implemented in this.
 ********************************************************************/
int main()
{
	int q,i;
	/******variables for Time evaluation for different version of algorithm**********************************/
	struct timeval ompStartTime,ompMidTime,ompEndTime;
	double ompGaussTime,ompBackstageTime,ompTotalTime;
	/*********************************************************************************************************

	/********Scan order of matrix from user*********/
	int n;
	printf("\nEnter the order of matrix: ");
	    scanf("%d",&n);
	/**********************************************/

	/*****Scan number of Thread from user******/
	int noThread;
	printf("\nEnter number of Thread: ");
	scanf("%d",&noThread);
	/*****************************************/


	/***********Dynamic Allowcation of matrix and X vector***********/
	double **A = (double **)calloc(n,sizeof(double*));
	for (q=0; q < n; q++)
		A[q] = (double*)calloc(n+1,sizeof(double*));

	/************Allocate Vector 'x for Pthread Answer8***************/
	double* x = (double*) malloc(sizeof(double)*n);
	/****************************************************************/

	t_info.n=n;
	t_info.numThreads=noThread;

    /*Generating random matrix by calling below function. */
	generateMatrix(n,A);


	t_info.A=A;
    /*printMarix("Initial Matrix :-",A);*/

/**************************************************************************************
 * OpenMP Algorith Implementation
 * Gaussian Elimination and Back Substitution and time calculation
 * *************************************************************************************/
	gettimeofday(&ompStartTime, NULL);
	openMPAlgorith(A);
	/*printMarix("After Matrix:-",A);*/
	gettimeofday(&ompMidTime, NULL);
	backSubstitution(A,x);
	gettimeofday(&ompEndTime, NULL);
/**************************************************************************************/


/******************Calculate Run Times and print out Solutions**********************************************/

		/***************OpenMP Time calculation and Print****************/
		ompGaussTime = ((ompMidTime.tv_sec  - ompStartTime.tv_sec) * 1000000u + ompMidTime.tv_usec - ompStartTime.tv_usec) / 1.e6;
		ompBackstageTime = ((ompEndTime.tv_sec  - ompMidTime.tv_sec) * 1000000u + ompEndTime.tv_usec - ompMidTime.tv_usec) / 1.e6;
		ompTotalTime = ((ompEndTime.tv_sec  - ompStartTime.tv_sec) * 1000000u + ompEndTime.tv_usec - ompStartTime.tv_usec) / 1.e6;

		printf("OpenMP Gaussian Elimination execution time: %.3f seconds.\n", ompGaussTime);
		printf("OpenMP Substitution execution time: %.3f seconds.\n", ompBackstageTime);
		printf("OpenMP Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", ompTotalTime,noThread);

		printf("\nThe solution is: \n");
		printf("X\t\tOpenMp\n");
		for(i=0; i<n; i++)
		{
			printf("x%d\t\t%.3f\n",i,x[i]); /* x1, x2, x3 are the required solutions*/
		}
	return(0);
} /* main() ends */


/***************************************************************
 * OpneMP algorithm.
 * This function is implemented for OPENMP Algorithm.
 **************************************************************/
void openMPAlgorith(double** A)
{
	int i,j,k;
	double c;
	int n=t_info.n;
	omp_set_nested(1);
	for(j=0; j<t_info.n-1; j++)  /*loop for the generation of upper triangular matrix*/
	{
		/*define parallel region*/
		#pragma omp parallel num_threads(t_info.numThreads) shared(n,A,j) private(i,c,k)
		{
			/*parallelizing for loop that has no data dependency.*/
			#pragma omp for schedule(static)
			for(i=j+1; i<n; i++)
			{
				c=A[i][j]/A[j][j];
				for(k=j; k<n+1; k++)
				{
					A[i][k]=A[i][k]-c*A[j][k];
				}
			}
		}
	}
}

