all: serialGaussianAlgorithm pThreadGaussianAlgorithm openMPGaussianAlgorithm pThreadPartialPivot pThreadGaussianVersion2 masterAllAlgrith  

sequential:
	
	gcc -o serialGaussianAlgorithm serialGaussianAlgorithm.c
pthread:
	
	gcc -o pThreadGaussianAlgorithm pThreadGaussianAlgorithm.c -lpthread
openmp:
	gcc -o openMPGaussianAlgorithm  openMPGaussianAlgorithm.c -lgomp -fopenmp 

pthreadPartialPivot:
	gcc -o pThreadPartialPivot pThreadPartialPivot.c -lpthread

pThreadGaussianVersion2:
	gcc -o pThreadGaussianVersion2 pThreadGaussianVersion2.c -lpthread

masterAllAlgrith:
	gcc -o masterAllAlgrith masterAllAlgrith.c -lpthread -lgomp -fopenmp

clean:
	rm serialGaussianAlgorithm	
	rm pThreadGaussianAlgorithm
	rm openMPGaussianAlgorithm
	rm pThreadPartialPivot
	rm pThreadGaussianVersion2
	rm masterAllAlgrith



