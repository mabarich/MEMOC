/**
 * @file main.cpp
 * @brief 
 */


#include <stdexcept>
#include <ctime>
#include <sys/time.h>

#include "TSPSolver.h"

// error status and messagge buffer
int status;
char errmsg[255];


int main (int argc, char const *argv[])
{
	try
	{
		if (argc < 3) throw std::runtime_error("usage: ./main #nodi disposizione euristica [tabuLength] [maxIter] [# clusters] solIniziale");
		TSP tspInstance;
		int cluster=1;
		int num=atoi(argv[1]);
		//0=4 punti, 1=random, 2=circonferenza, 3=clusters, 4=leggi da file
		int fun=atoi(argv[2]); 
		//0=LS, 1=TS, 2=ATS
		int met=atoi(argv[3]);
		int tabuLength = 8;                                                           
		int maxIter    = 110;
		if (met==1 || met==2)
		{
			tabuLength = atoi(argv[4]);                                                           
			maxIter    = atoi(argv[5]);
		}
		//Secondo argomento è il numero di cluster
		if (fun==3)
			cluster=atoi(argv[6]);
		tspInstance.read(num, fun, cluster);
		TSPSolution aSolution(tspInstance);
		//Tempo di inizio
		clock_t t1,t2;
		t1 = clock();
		struct timeval  tv1, tv2;
		gettimeofday(&tv1, NULL);
		//Soluzione iniziale random
		TSPSolver tspSolver;
		int iniz=atoi(argv[7]);
		if (iniz==0)
			tspSolver.initRnd(aSolution);
		else if (iniz==1)
			tspSolver.initBst(aSolution, tspInstance);
		//Ricerca dei vicinati
		TSPSolution bestSolution(tspInstance);
		if (met==0)
			tspSolver.solve(tspInstance,aSolution,bestSolution);
		else if (met==1)
			tspSolver.solveTL(tspInstance,aSolution,tabuLength,maxIter,bestSolution);
		else if (met==2)
			tspSolver.solveAS(tspInstance,aSolution,tabuLength,maxIter,bestSolution);
		//Tempo di fine
		t2 = clock();
		gettimeofday(&tv2, NULL);
		//Stampe
		std::cout << "FROM solution: "; 
		aSolution.print();
		std::cout << "(value : " << tspSolver.evaluate(aSolution,tspInstance) << ")\n";
		std::cout << "TO   solution: "; 
		bestSolution.print();
		std::cout << "(value : " << tspSolver.evaluate(bestSolution,tspInstance) << ")\n";
		std::cout << "in " << (double)(tv2.tv_sec+tv2.tv_usec*1e-6 - (tv1.tv_sec+tv1.tv_usec*1e-6)) << " seconds (user time)\n";
		std::cout << "in " << (double)(t2-t1) / CLOCKS_PER_SEC << " seconds (CPU time)\n";
		
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}