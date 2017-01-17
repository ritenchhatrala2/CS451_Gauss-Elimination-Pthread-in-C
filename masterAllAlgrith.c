/*
 * pgauss.c
 *
 *  Created on: Sep 28, 2016
 *      Author: riten
 */
#include <stdio.h>      /* Input/Output */
#include <stdlib.h>     /* General Utilities */
#include <pthread.h>    /* POSIX Threads */
#include <omp.h>		/* OpenMP Lib*/

/*float A[4][5]={
		{2,1,-1,2,5},
		{4,5,-3,6,9},
		{-2,5,-2,6,4},
		{4,11,-4,8,2}};*///ans 3 1 -2 1

/*float A[3][4]={
		{3,2,-4,3},
		{2,3,3,15},
		{5,-3,1,14}};*/

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
void pfunction ( void *ptr );
void pfunctionV2 ( void *ptr );
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

/****************************************************************
 * Serial Gaussian Elimination Algorithm implementation
 ****************************************************************/
void serialAlgorith(double** A){
	int i,j,k;

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

void partialPivot(double** A){

	int   i,k,v;
	int n=t_info.n;
	int j=t_info.j;
	double checkmax, temp, max;
	max = (double)fabs(A[j][j]);
    v = j;
    /* Searching for row with largest pivot */
    for (i=j+1; i<n; i++){
    	checkmax = (double) fabs(A[i][j]);
    	if(checkmax > max) {max = checkmax; v=i;}
    }
    /* Row interchanges */
    if(v != j) {
    	for(k=j; k<n+1; k++) {
    		temp = A[j][k];
    		A[j][k] = A[v][k];
    		A[v][k] = temp;
    	}
    }
    /*printMarix("No j %d:-",A);*/
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
	struct timeval serialStartTime,serialEndTime,serialTotalTime,pThreadStartTime, pThreadMidTime, pThreadEndTime;
	double pthreadGaussTime, pthreadBackStageTime, pthreadTotalTime,serialGausstime,serialBackstageTime,serialTTime;

	struct timeval ompStartTime,ompMidTime,ompEndTime;
	double ompGaussTime,ompBackstageTime,ompTotalTime;

	struct timeval pThreadStartTimeV2, pThreadMidTimeV2, pThreadEndTimeV2;
	double pthreadGaussTimeV2, pthreadBackStageTimeV2, pthreadTotalTimeV2;

	struct timeval pThreadStartTimePP, pThreadMidTimePP, pThreadEndTimePP;
	double pthreadGaussTimePP, pthreadBackStageTimePP, pthreadTotalTimePP;
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

	pthread_t thread[noThread];

	/***********Dynamic Allowcation of matrix and X vector***********/
	double **A = (double **)calloc(n,sizeof(double*));
	for (q=0; q < n; q++)
		A[q] = (double*)calloc(n+1,sizeof(double*));

	/************Allocate Vector 'x for Pthread Answer8***************/
	double* x = (double*) malloc(sizeof(double)*n);
	/****************************************************************/

	/*Allocate Vector 'y for Serial Algorithm)*/
	double* y = (double*) malloc(sizeof(double)*n);
	/****************************************************************/


	/*Allocate Vector 'z for OpenMP algorithm.*/
	double* z = (double*) malloc(sizeof(double)*n);
	/****************************************************************/

	/*Allocate Vector 'w for pthreadV2 algorithm.*/
	double* w = (double*) malloc(sizeof(double)*n);
	/****************************************************************/

	/*Allocate Vector 'pp for pthread partial pivot algorithm.*/
	double* pp = (double*) malloc(sizeof(double)*n);
	/****************************************************************/

    t_info.n=n;
    t_info.numThreads=noThread;

    /*Generating random matrix by calling below function. */
	generateMatrix(n,A);


	t_info.A=A;
	double** B=A;
	double** C=A;
	double** P=A;
    /*printMarix("Initial Matrix :-",A);*/

	/*******initialization of thread_index array to pass thread index in pthread fuction******/
	int *thread_index = calloc (n, sizeof (int));
	for(i = 0; i < noThread; i++)
	{
		thread_index[i] = i;
	}
	/*****************************************************************************************/

/**************************************************************************************
 * Pthread Algorithm.
 * Gaussian Elimination and Back Substitution for PThread and time evaluation for PThread
 * *************************************************************************************/
	gettimeofday(&pThreadStartTime, NULL);

	/*Gaussian Elimination*/
	for(j=0; j<n-1; j++)
	{
		t_info.j=j;
		for(d=0;d<noThread;d++)
		{
			pthread_create (&thread[d], NULL, (void *) &pfunction, (void*)&thread_index[d]);
		}
		for(l=0; l<noThread; l++)
		{
			pthread_join(thread[l], NULL);
		}
	}

	gettimeofday(&pThreadMidTime, NULL);

	/*printMarix("After Matrix:-",A);*/

	/*back Substitution*/
	backSubstitution(A,x);

	gettimeofday(&pThreadEndTime, NULL);

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
	backSubstitution(A,w);

	gettimeofday(&pThreadEndTimeV2, NULL);

/**************************************************************************************
 * Serial Algorith Implementation
 * Gaussian Elimination and Back Substitution and time calculation
 * *************************************************************************************/

	gettimeofday(&serialStartTime, NULL);
	/*Gaussian Elimination*/
	serialAlgorith(B);
	gettimeofday(&serialEndTime, NULL);
	/*back Substitution*/
	backSubstitution(B,y);
	gettimeofday(&serialTotalTime, NULL);

/**************************************************************************************
 * OpenMP Algorith Implementation
 * Gaussian Elimination and Back Substitution and time calculation
 * *************************************************************************************/
	gettimeofday(&ompStartTime, NULL);
	openMPAlgorith(C);
	gettimeofday(&ompMidTime, NULL);
	backSubstitution(C,z);
	gettimeofday(&ompEndTime, NULL);
/**************************************************************************************/

/**************************************************************************************
 * Pthread with partial pivot Algorithm.
 * Gaussian Elimination and Back Substitution for PThread and time evaluation for PThread
 * *************************************************************************************/
	gettimeofday(&pThreadStartTimePP, NULL);

	/*Gaussian Elimination*/
	for(j=0; j<n-1; j++) /* loop for the generation of upper triangular matrix*/
	{

		t_info.j=j;
		partialPivot(P);
		for(d=0;d<noThread;d++)
		{
			pthread_create (&thread[d], NULL, (void *) &pfunction, (void*)&thread_index[d]);
		}
		for(l=0; l<noThread; l++) /* loop for the generation of upper triangular matrix*/
		{
			pthread_join(thread[l], NULL);
		}
	}

	gettimeofday(&pThreadMidTimePP, NULL);

	/*printMarix("After Matrix:-",A);*/

	/*back Substitution*/
	backSubstitution(P,pp);

	gettimeofday(&pThreadEndTimePP, NULL);



/******************Calculate Run Times and print out Solutions**********************************************/

		/***************Serial Algorithm Time calculation and Print****************/
		serialGausstime = ((serialEndTime.tv_sec  - serialStartTime.tv_sec) * 1000000u + serialEndTime.tv_usec - serialStartTime.tv_usec) / 1.e6;
		serialBackstageTime = ((serialTotalTime.tv_sec  - serialEndTime.tv_sec) * 1000000u + serialTotalTime.tv_usec - serialEndTime.tv_usec) / 1.e6;
		serialTTime = ((serialTotalTime.tv_sec  - serialStartTime.tv_sec) * 1000000u + serialTotalTime.tv_usec - serialStartTime.tv_usec) / 1.e6;

		printf("Serial Gaussian Elimination execution time: %.3f seconds.\n", serialGausstime);
		printf("Serial Back Substitution execution time: %.3f seconds.\n", serialBackstageTime);
		printf("Total serial execution: \n%.3f seconds elapsed.\n\n", serialTTime);

		/***************Pthread Time calculation and Print****************/
		pthreadGaussTime = ((pThreadMidTime.tv_sec  - pThreadStartTime.tv_sec) * 1000000u + pThreadMidTime.tv_usec - pThreadStartTime.tv_usec) / 1.e6;
		pthreadBackStageTime = ((pThreadEndTime.tv_sec  - pThreadMidTime.tv_sec) * 1000000u + pThreadEndTime.tv_usec - pThreadMidTime.tv_usec) / 1.e6;
		pthreadTotalTime = ((pThreadEndTime.tv_sec  - pThreadStartTime.tv_sec) * 1000000u + pThreadEndTime.tv_usec - pThreadStartTime.tv_usec) / 1.e6;

		printf("Pthread Gaussian Elimination execution time: %.3f seconds.\n", pthreadGaussTime);
		printf("Pthread Back Substitution execution time: %.3f seconds.\n", pthreadBackStageTime);
		printf("Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", pthreadTotalTime, noThread);

		/***************OpenMP Time calculation and Print****************/
		ompGaussTime = ((ompMidTime.tv_sec  - ompStartTime.tv_sec) * 1000000u + ompMidTime.tv_usec - ompStartTime.tv_usec) / 1.e6;
		ompBackstageTime = ((ompEndTime.tv_sec  - ompMidTime.tv_sec) * 1000000u + ompEndTime.tv_usec - ompMidTime.tv_usec) / 1.e6;
		ompTotalTime = ((ompEndTime.tv_sec  - ompStartTime.tv_sec) * 1000000u + ompEndTime.tv_usec - ompStartTime.tv_usec) / 1.e6;

		printf("OpenMP Gaussian Elimination execution time: %.3f seconds.\n", ompGaussTime);
		printf("OpenMP Substitution execution time: %.3f seconds.\n", ompBackstageTime);
		printf("OpenMP Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", ompTotalTime, noThread);


		/***************Version-2 Pthread Time calculation and Print***************/
		pthreadGaussTimeV2 = ((pThreadMidTimeV2.tv_sec  - pThreadStartTimeV2.tv_sec) * 1000000u + pThreadMidTimeV2.tv_usec - pThreadStartTimeV2.tv_usec) / 1.e6;
		pthreadBackStageTimeV2 = ((pThreadEndTimeV2.tv_sec  - pThreadMidTimeV2.tv_sec) * 1000000u + pThreadEndTimeV2.tv_usec - pThreadMidTimeV2.tv_usec) / 1.e6;
		pthreadTotalTimeV2 = ((pThreadEndTimeV2.tv_sec  - pThreadStartTimeV2.tv_sec) * 1000000u + pThreadEndTimeV2.tv_usec - pThreadStartTimeV2.tv_usec) / 1.e6;

		printf("Version2 Pthread Gaussian Elimination execution time: %.3f seconds.\n", pthreadGaussTimeV2);
		printf("Version2 Pthread Back Substitution execution time: %.3f seconds.\n", pthreadBackStageTimeV2);
		printf("Version2 Pthread Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", pthreadTotalTimeV2, noThread);


		/***************Partial Pivot Pthread Time calculation and Print***************/
		pthreadGaussTimePP = ((pThreadMidTimePP.tv_sec  - pThreadStartTimePP.tv_sec) * 1000000u + pThreadMidTimePP.tv_usec - pThreadStartTimePP.tv_usec) / 1.e6;
		pthreadBackStageTimePP = ((pThreadEndTimePP.tv_sec  - pThreadMidTimePP.tv_sec) * 1000000u + pThreadEndTimePP.tv_usec - pThreadMidTimePP.tv_usec) / 1.e6;
		pthreadTotalTimePP = ((pThreadEndTimePP.tv_sec  - pThreadStartTimePP.tv_sec) * 1000000u + pThreadEndTimePP.tv_usec - pThreadStartTimePP.tv_usec) / 1.e6;

		printf("Partial Pivot Pthread Gaussian Elimination execution time: %.3f seconds.\n", pthreadGaussTimePP);
		printf("Partial Pivot Pthread Back Substitution execution time: %.3f seconds.\n", pthreadBackStageTimePP);
		printf("Partial Pivot Pthread Total execution: \n%.3f seconds elapsed with %d threads used.\n\n", pthreadTotalTimePP, noThread);


		printf("\nTotal Time for all Algorithm:-\n");
		printf("_____________________________________________________________________________________________________\n");
		printf("|core\t\t|serial\t\t|pThread\t|pThreadV2\t|OpenMP\t\t|PartialPivotPthread|\n");
		printf("_____________________________________________________________________________________________________\n");
		printf("|%d\t\t|%.3f\t\t|%.3f\t\t|%.3f\t\t|%.3f\t\t|%.3f\t\t    |\n",noThread,serialTTime,pthreadTotalTime,
				pthreadTotalTimeV2,ompTotalTime,pthreadTotalTimePP);
		printf("-----------------------------------------------------------------------------------------------------\n");

		printf("\nThe solution is: \n");
		printf("X\t\tSerial\t\tPthread\t\tOpenMp\t\tpthreadV2\tpartialPivotPthread\n");
		for(i=0; i<n; i++)
		{
			printf("x%d\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n",i,y[i],x[i],z[i],w[i],pp[i]); /* x1, x2, x3 are the required solutions*/
		}
	return(0);
} /* main() ends */


/***************************************************************
 * PThread Parallel Function.
 * This function is executed by every thread creaded in main().
 **************************************************************/
void pfunction ( void *ptr )
{
	int i,k;
    int thread_index= *(int *)ptr;
    int j=t_info.j;
    double c;
    double** A=t_info.A;

    for(i=j+1+thread_index; i<t_info.n; i=i+t_info.numThreads)
	{
		c=A[i][j]/A[j][j];
		for(k=j; k<t_info.n+1; k++)
		{
			A[i][k]=A[i][k]-c*A[j][k];
		}
	}
    pthread_exit(0); /* exit */
}

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
			A[i][k]=A[i][k]-c*A[j][k];
		}

    pthread_exit(0); /* exit */
}

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
	for(j=0; j<t_info.n-1; j++)
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
