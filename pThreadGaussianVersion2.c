/*
 * pthreadV2.c
 *
 *  Created on: Oct 2, 2016
 *      Author: riten
 */
/*
 * pthread.c
 *
 *  Created on: Oct 1, 2016
 *      Author: riten
 */
#include <stdio.h>      /* Input/Output */
#include <stdlib.h>     /* General Utilities */
#include <pthread.h>    /* POSIX Threads */

/*************************************************
 * Structure to globalize the data for Thread Use.
 *************************************************/
struct Thread_Data {
	double** A;
	int n;
	int j;
	int i;
	double c;
	int numThreads;
} t_info;

/*Declareration of Pthread and openMP function*/
void pfunctionV2 ( void *ptr );

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
	/*A[0][0]=3;
		A[0][1]=2;
		A[0][2]=-4;
		A[0][3]=3;
		A[1][0]=2;
		A[1][1]=3;
		A[1][2]=3;
		A[1][3]=15;
		A[2][0]=5;
		A[2][1]=-3;
		A[2][2]=1;
		A[2][3]=14;*/
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

/*******************************************************************
 * Main() function execution start.
 * Serial, Pthread , Pthread Version 2 and openMP implemented in this.
 ********************************************************************/
int main()
{
	int q,i,j,d,l;
	double c;
	/******variables for Time evaluation for different version of algorithm**********************************/
	struct timeval pThreadStartTimeV2, pThreadMidTimeV2, pThreadEndTimeV2;
	double pthreadGaussTimeV2, pthreadBackStageTimeV2, pthreadTotalTimeV2;

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

	pthread_t thread[noThread];

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
    printMarix("Initial Matrix :-",A);

	/*******initialization of thread_index array to pass thread index in pthread fuction******/
	int *thread_index = calloc (n, sizeof (int));
	for(i = 0; i < noThread; i++)
	{
		thread_index[i] = i;
	}
	/*****************************************************************************************/

	/**************************************************************************************
	 * Pthread Algorithm Version 2.
	 * Gaussian Elimination and Back Substitution for PThread and time evaluation for PThread
	 * *************************************************************************************/
		gettimeofday(&pThreadStartTimeV2, NULL);

		/*Gaussian Elimination*/
		for(j=0; j<n-1; j++)
		{
			t_info.j=j;
			for(i=j+1; i<t_info.n; i++)
			{
				c=A[i][j]/A[j][j];
				t_info.i=i;
				t_info.c=c;
				for(d=0;d<noThread;d++)
				{
					pthread_create (&thread[d], NULL, (void *) &pfunctionV2, (void*)&thread_index[d]);
				}
				for(l=0; l<noThread; l++)
				{
					pthread_join(thread[l], NULL);
				}
			}
		}

		gettimeofday(&pThreadMidTimeV2, NULL);

		/*printMarix("After Matrix:-",A);*/

		/*back Substitution*/
		backSubstitution(A,x);

		gettimeofday(&pThreadEndTimeV2, NULL);

/******************Calculate Run Times and print out Solutions**********************************************/

		/***************Version-2 Pthread Time calculation and Print***************/
		pthreadGaussTimeV2 = ((pThreadMidTimeV2.tv_sec  - pThreadStartTimeV2.tv_sec) * 1000000u + pThreadMidTimeV2.tv_usec - pThreadStartTimeV2.tv_usec) / 1.e6;
		pthreadBackStageTimeV2 = ((pThreadEndTimeV2.tv_sec  - pThreadMidTimeV2.tv_sec) * 1000000u + pThreadEndTimeV2.tv_usec - pThreadMidTimeV2.tv_usec) / 1.e6;
		pthreadTotalTimeV2 = ((pThreadEndTimeV2.tv_sec  - pThreadStartTimeV2.tv_sec) * 1000000u + pThreadEndTimeV2.tv_usec - pThreadStartTimeV2.tv_usec) / 1.e6;

		printf("Version2 Pthread Gaussian Elimination execution time: %.3f seconds.\n", pthreadGaussTimeV2);
		printf("Version2 Pthread Back Substitution execution time: %.3f seconds.\n", pthreadBackStageTimeV2);
		printf("Version2 Pthread Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", pthreadTotalTimeV2, noThread);

		printf("\nThe solution is: \n");
		printf("X\t\tPthread\n");
		for(i=0; i<n; i++)
		{
			printf("x%d\t\t%.3f\n",i,x[i]); /* x1, x2, x3 are the required solutions*/
		}
	return(0);
} /* main() ends */

/***************************************************************
 * PThread version 2 Parallel Function.
 * This function is executed by every thread creaded in main().
 **************************************************************/
void pfunctionV2 ( void *ptr )
{
	int k;
    int thread_index= *(int *)ptr;
    int j=t_info.j;
    int i=t_info.i;
    double c=t_info.c;
    double** A=t_info.A;

		for(k=j+thread_index; k<t_info.n+1; k=k+t_info.numThreads)
		{
			printf("\nj=%d\ti=%d\tc=%d",j,i,c);
			A[i][k]=A[i][k]-c*A[j][k];
		}

    pthread_exit(0); /* exit */
}
