/*
 * serial.c
 *
 *  Created on: Oct 1, 2016
 *      Author: riten
 */
#include <stdio.h>      /* Input/Output */
#include <stdlib.h>     /* General Utilities */

/*************************************************
 * Structure to globalize the data for Thread Use.
 *************************************************/
struct Thread_Data {
	double** A;
	int n;
	int numThreads;
} t_info;


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
				printf("%lf\t",A[i][j]);
			}
			printf("\n");
		}
}

/************************************************************
 * Back Substitution implementation after matrix transform to
 * upper triangular form
 ************************************************************/
void backSubstitution(double** A,double* x){
	double sum;
	int i,j;
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

/****************************************************************
 * Serial Gaussian Elimination Algorithm implementation
 ****************************************************************/
void serialAlgorith(double** A){

	int j,i,k;

	for(j=0; j<t_info.n-1; j++) /* loop for the generation of upper triangular matrix*/
	{
		for(i=j+1; i<t_info.n; i++)
		{
			double c=A[i][j]/A[j][j];
			for(k=j; k<t_info.n+1; k++)
			{
				A[i][k]=A[i][k]-c*A[j][k];
			}
		}
	}
}

/*******************************************************************
 * Main() function execution start.
 * Serial implemented in this.
 ********************************************************************/
int main()
{
	int q,i;
	/******variables for Time evaluation for different version of algorithm**********************************/
	struct timeval serialStartTime,serialEndTime,serialTotalTime;
	double serialGausstime,serialBackstageTime,serialTTime;
	/*********************************************************************************************************

	/********Scan order of matrix from user*********/
	int n;
	printf("\nEnter the order of matrix: ");
	    scanf("%d",&n);
	/**********************************************/

	/***********Dynamic Allowcation of matrix and X vector***********/
	double **A = (double **)calloc(n,sizeof(double*));
	for (q=0; q < n; q++)
		A[q] = (double*)calloc(n+1,sizeof(double*));

	/************Allocate Vector 'y for Serial Algorithm)***************/
	double* x = (double*) malloc(sizeof(double)*n);
	/****************************************************************/

    t_info.n=n;

    /*Generating random matrix by calling below function. */
	generateMatrix(n,A);


	t_info.A=A;
    /*printMarix("Initial Matrix :-",A);*/

/**************************************************************************************
 * Serial Algorith Implementation
 * Gaussian Elimination and Back Substitution and time calculation
 * *************************************************************************************/

	gettimeofday(&serialStartTime, NULL);
	/*Gaussian Elimination*/
	serialAlgorith(A);
	gettimeofday(&serialEndTime, NULL);
	/*back Substitution*/
	backSubstitution(A,x);
	gettimeofday(&serialTotalTime, NULL);

/******************Calculate Run Times and print out Solutions**********************************************/

		/***************Serial Algorithm Time calculation and Print****************/
		serialGausstime = ((serialEndTime.tv_sec  - serialStartTime.tv_sec) * 1000000u + serialEndTime.tv_usec - serialStartTime.tv_usec) / 1.e6;
		serialBackstageTime = ((serialTotalTime.tv_sec  - serialEndTime.tv_sec) * 1000000u + serialTotalTime.tv_usec - serialEndTime.tv_usec) / 1.e6;
		serialTTime = ((serialTotalTime.tv_sec  - serialStartTime.tv_sec) * 1000000u + serialTotalTime.tv_usec - serialStartTime.tv_usec) / 1.e6;

		printf("Serial Gaussian Elimination execution time: %.3f seconds.\n", serialGausstime);
		printf("Serial Back Substitution execution time: %.3f seconds.\n", serialBackstageTime);
		printf("Total serial execution: \n%.3f seconds elapsed.\n\n", serialTTime);

		printf("\nThe solution is: \n");
		printf("X\t\tSerial\n");
		for(i=0; i<n; i++)
		{
			printf("x%d\t\t%.3f\n",i,x[i]); /* x1, x2, x3 are the required solutions*/
		}
	return(0);
} /* main() ends */

