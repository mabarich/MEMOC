/**
 * @file TSPSolver.cpp
 * @brief TSP solver (neighborhood search)
 *
 */

#include "TSPSolver.h"
#include <iostream>


bool TSPSolver::solve ( const TSP& tsp , const TSPSolution& initSol , TSPSolution& bestSol, int opt )
{
	try
	{
		bool stop = false;
		int  iter = 0;
		//Valuto la soluzione corrente
		TSPSolution currSol(initSol);
		double bestValue, currValue;
		bestValue = currValue = evaluate(currSol,tsp);
		TSPMove move;
		while ( ! stop) 
		{
			if ( tsp.n < 20 ) currSol.print();
			std::cout << " (" << ++iter << ") value " << currValue << " (" << evaluate(currSol,tsp) << ")";
			//Passa ogni vicino e trova il miglior risultato
			double bestNeighValue;
			if (opt==0)
			{	
				bestNeighValue= currValue + findBestNeighbor(tsp,currSol,move);
				std::cout << " move: " << move.from << " , " << move.to << std::endl;
			}
			else if (opt==1)
			{
				bestNeighValue= currValue+ findBestNeighbor3opt(tsp,currSol,move);
				std::cout << " move: " << move.from << " , " << move.to <<  " , " << move.d << std::endl;
			}
			//Criterio di stop: se ho un risultato migliorante, vado la, altrimenti mi fermo
			if ( bestNeighValue < currValue )
			{
				bestValue = currValue = bestNeighValue;
				//Passo al vicino
				if (opt==0)
					currSol = swap(currSol,move);
				else if (opt==1)
					currSol = swap3opt(currSol,move);
				stop = false;
			} 
			else 
			{
				cout<<"Non trovo di meglio";
				stop = true;
			}
		}
		bestSol = currSol; 
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
		return false;
	}
	return true;
}

bool TSPSolver::solveTL ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol, int opt )
{
	try
	{
		bool stop = false;
		int  iter = 0;
		//Tabu List
		tabuLength = tabulength;
		tabuList.reserve(tsp.n);
		initTabuList(tsp.n);
		//Valuto la soluzione corrente
		TSPSolution currSol(initSol);
		double bestValue, currValue;
		bestValue = currValue = evaluate(currSol,tsp);
		TSPMove move;
		while ( ! stop ) 
		{
			if ( tsp.n < 20 ) currSol.print();
			std::cout << " (" << ++iter << ") value " << currValue << " (" << evaluate(currSol,tsp) << ")";
			int bestNeighValue;
			//Passa ogni vicino e trova il miglior risultato
			if (opt==0)
			{	
				bestNeighValue= currValue + findBestNeighbor2(tsp,currSol,iter,move);
				std::cout << " move: " << move.from << " , " << move.to << std::endl;
			}
			else if (opt==1)
			{
				bestNeighValue= currValue+ findBestNeighbor3opt2(tsp,currSol,iter,move);
				std::cout << " move: " << move.from << " , " << move.to <<  " , " << move.d << std::endl;
			}
			//Criterio di stop: se non ho più vicini legali mi fermo
			if ( bestNeighValue >= tsp.infinite ) 
			{                                                             
				std::cout << "\tmove: NO legal neighbour" << std::endl;                                           
				stop = true;                                                                                      
				continue;                                                                                         
			}
			
			std::cout << "\tmove: " << move.from << " , " << move.to;
			//Nella tabu list metto il numero dell'iterazione in cui è stata eseguita una mossa. Se la differenza tra l'iterazione corrente e quella salvata è maggiore della lunghezza della lista allora posso rifarla.
			tabuList[currSol.sequence[move.from]] = iter;                                                       
			tabuList[currSol.sequence[move.to]]   = iter;                                                       
			if (opt==0)
				currSol = swap(currSol,move);
			else if (opt==1)
				currSol = swap3opt(currSol,move);                                                                       
			currValue = bestNeighValue;                                                                          
			if ( currValue < bestValue - 0.01 ) 
			{                                                               
				bestValue = currValue;                                                                            
				bestSol = currSol;                                                                                
				std::cout << "\t***";                                                                             
			}                                                                                                   
			//Criterio di stop: se ho superato il massimo numero di iterazioni consentite mi fermo
			if ( iter > maxIter ) 
			{                                                                            
				stop = true;                                                                                      
			}                                                                                                   
			std::cout << std::endl;
			
		}
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
		return false;
	}
	return true;
}

/*bool TSPSolver::solveTL ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol, int opt )
{
	try
	{
		bool stop = false;
		int  iter = 0;
		//Tabu List
		tabuLength = tabulength;
		tabuList.reserve(tsp.n);
		initTabuList(tsp.n);
		//Valuto la soluzione corrente
		TSPSolution currSol(initSol);
		double bestValue, currValue;
		bestValue = currValue = evaluate(currSol,tsp);
		TSPMove move;
		while ( ! stop ) 
		{
			if ( tsp.n < 20 ) currSol.print();
			std::cout << " (" << ++iter << ") value " << currValue << " (" << evaluate(currSol,tsp) << ")";
			//Passa ogni vicino e trova il miglior risultato
			double bestNeighValue = currValue + findBestNeighbor2(tsp,currSol,iter,move);
			std::cout << " move: " << move.from << " , " << move.to << std::endl;
			//Criterio di stop: se non ho più vicini legali mi fermo
			if ( bestNeighValue >= tsp.infinite ) 
			{                                                             
				std::cout << "\tmove: NO legal neighbour" << std::endl;                                           
				stop = true;                                                                                      
				continue;                                                                                         
			}
			
			std::cout << "\tmove: " << move.from << " , " << move.to;
			//Nella tabu list metto il numero dell'iterazione in cui è stata eseguita una mossa. Se la differenza tra l'iterazione corrente e quella salvata è maggiore della lunghezza della lista allora posso rifarla.
			tabuList[currSol.sequence[move.from]] = iter;                                                       
			tabuList[currSol.sequence[move.to]]   = iter;                                                       
			currSol = swap(currSol,move);                                                                       
			currValue = bestNeighValue;                                                                          
			if ( currValue < bestValue - 0.01 ) 
			{                                                               
				bestValue = currValue;                                                                            
				bestSol = currSol;                                                                                
				std::cout << "\t***";                                                                             
			}                                                                                                   
			//Criterio di stop: se ho superato il massimo numero di iterazioni consentite mi fermo
			if ( iter > maxIter ) 
			{                                                                            
				stop = true;                                                                                      
			}                                                                                                   
			std::cout << std::endl;
			
		}
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
		return false;
	}
	return true;
}*/

/*bool TSPSolver::solveAS ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol, int opt )   
{
	try
	{
		bool stop = false;
		int  iter = 0;

		///Tabu Search
		tabuLength = tabulength;
		tabuList.reserve(tsp.n);
		initTabuList(tsp.n);
		///

		TSPSolution currSol(initSol);
		double bestValue, currValue;
		bestValue = currValue = evaluate(currSol,tsp);
		TSPMove move;
		while ( ! stop ) 
		{
			++iter;                                                                                             
			if ( tsp.n < 20 ) currSol.print();
			std::cout << " (" << iter << ") value " << currValue << "\t(" << evaluate(currSol,tsp) << ")";      
			double aspiration = bestValue-currValue;                                                            
			double bestNeighValue = currValue + findBestNeighbor3(tsp,currSol,iter,aspiration,move);             

			if ( bestNeighValue >= tsp.infinite ) 
			{                                                             
				std::cout << "\tmove: NO legal neighbour" << std::endl;                                           
				stop = true;                                                                                      
				continue;                                                                                         
			}                                                                                                   

			std::cout << "\tmove: " << move.from << " , " << move.to;

			tabuList[currSol.sequence[move.from]] = iter;                                                       
			tabuList[currSol.sequence[move.to]]   = iter;                                                       
			currSol = swap(currSol,move);                                                                       
			currValue = bestNeighValue;                                                                          
			if ( currValue < bestValue -0.01 ) 
			{                                                                
				bestValue = currValue;                                                                            
				bestSol = currSol;                                                                                
				std::cout << "\t***";                                                                             
			}                                                                                                   

			if ( iter > maxIter ) 
			{                                                                             
				stop = true;                                                                                      
			}                                                                                                  
			std::cout << std::endl;  
		}                                                                                 
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
		return false;
	}
	return true;
}*/

bool TSPSolver::solveAS ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol, int opt )   
{
	try
	{
		bool stop = false;
		int  iter = 0;

		///Tabu Search
		tabuLength = tabulength;
		tabuList.reserve(tsp.n);
		initTabuList(tsp.n);
		///

		TSPSolution currSol(initSol);
		double bestValue, currValue;
		bestValue = currValue = evaluate(currSol,tsp);
		TSPMove move;
		while ( ! stop ) 
		{
			++iter;                                                                                             
			if ( tsp.n < 20 ) currSol.print();
			std::cout << " (" << iter << ") value " << currValue << "\t(" << evaluate(currSol,tsp) << ")";      
			double aspiration = bestValue-currValue;                                                             
			double bestNeighValue;
			//Passa ogni vicino e trova il miglior risultato
			if (opt==0)
			{	
				bestNeighValue= currValue + findBestNeighbor2(tsp,currSol,iter,move);
				std::cout << " move: " << move.from << " , " << move.to << std::endl;
			}
			else if (opt==1)
			{
				bestNeighValue= currValue+ findBestNeighbor3opt3(tsp,currSol,iter,aspiration,move);
				std::cout << " move: " << move.from << " , " << move.to <<  " , " << move.d << std::endl;
			}
			         

			if ( bestNeighValue >= tsp.infinite ) 
			{                                                             
				std::cout << "\tmove: NO legal neighbour" << std::endl;                                           
				stop = true;                                                                                      
				continue;                                                                                         
			}                                                                                                   

			std::cout << "\tmove: " << move.from << " , " << move.to;

			tabuList[currSol.sequence[move.from]] = iter;                                                       
			tabuList[currSol.sequence[move.to]]   = iter;                                                       
			if (opt==0)
				currSol = swap(currSol,move);
			else if (opt==1)
				currSol = swap3opt(currSol,move);                                                                       
			currValue = bestNeighValue;                                                                          
			if ( currValue < bestValue -0.01 ) 
			{                                                                
				bestValue = currValue;                                                                            
				bestSol = currSol;                                                                                
				std::cout << "\t***";                                                                             
			}                                                                                                   

			if ( iter > maxIter ) 
			{                                                                             
				stop = true;                                                                                      
			}                                                                                                  
			std::cout << std::endl;  
		}                                                                                 
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
		return false;
	}
	return true;
}

//2-opt. Scambia due elementi non consecutivi della sequenza
TSPSolution& TSPSolver::swap ( TSPSolution& tspSol , const TSPMove& move ) 
{
	TSPSolution tmpSol(tspSol);
	for ( int i = move.from ; i <= move.to ; ++i ) 
	{
		tspSol.sequence[i] = tmpSol.sequence[move.to-(i-move.from)];
	}
	return tspSol;
}

//Trova la minor differenza, non il vicino
double TSPSolver::findBestNeighbor ( const TSP& tsp , const TSPSolution& currSol , TSPMove& move )
/* Determine the *move* yielding the best 2-opt neigbor solution
 */
{
	double bestDecrement = tsp.infinite;
	// intial and final position are fixed (initial/final node remains 0)
	for ( uint a = 1 ; a < currSol.sequence.size() - 2 ; a++ ) 
	{
		int h = currSol.sequence[a-1];
		int i = currSol.sequence[a];
		for ( uint b = a + 1 ; b < currSol.sequence.size() - 1 ; b++ ) 
		{
			int j = currSol.sequence[b];
			int l = currSol.sequence[b+1];
			double neighDecrement = - tsp.cost[h][i] - tsp.cost[j][l] + tsp.cost[h][j] + tsp.cost[i][l] ;
			if ( neighDecrement < bestDecrement ) 
			{
				bestDecrement = neighDecrement;
				move.from = a;
				move.to = b;
			}
		}
	}
	  
	return bestDecrement;
}

//Trova la minor differenza, non il vicino con la tabu list
double TSPSolver::findBestNeighbor2 ( const TSP& tsp , const TSPSolution& currSol , int currIter , TSPMove& move )
/* Determine the *move* yielding the best 2-opt neigbor solution
 */
{
	double bestDecrement = tsp.infinite;
	// intial and final position are fixed (initial/final node remains 0)
	for ( uint a = 1 ; a < currSol.sequence.size() - 2 ; a++ ) 
	{
		int h = currSol.sequence[a-1];
		int i = currSol.sequence[a];
		for ( uint b = a + 1 ; b < currSol.sequence.size() - 1 ; b++ ) 
		{
			int j = currSol.sequence[b];
			int l = currSol.sequence[b+1];	
			if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) ) continue;
			double neighDecrement = - tsp.cost[h][i] - tsp.cost[j][l] + tsp.cost[h][j] + tsp.cost[i][l] ;
			if ( neighDecrement < bestDecrement ) 
			{
				bestDecrement = neighDecrement;
				move.from = a;
				move.to = b;
			}
		}
	}
	  
	return bestDecrement;
}


double TSPSolver::findBestNeighbor3 ( const TSP& tsp , const TSPSolution& currSol , int currIter , double aspiration , TSPMove& move )     //**// TSAC: use aspiration
/* Determine the NON-TABU (or satisfying aspiration) *move* yielding the best 2-opt neigbor solution
 * Aspiration criteria: 'neighDecrement' better than 'aspiration' (notice that 'aspiration'
 * has been set such that if 'neighDecrement' is better than 'aspiration' than we have a
 * nee incumbent solution)
 */
{
	double bestDecrement = tsp.infinite;

	for ( uint a = 1 ; a < currSol.sequence.size() - 2 ; a++ ) 
	{
		int h = currSol.sequence[a-1];
		int i = currSol.sequence[a];
		for ( uint b = a + 1 ; b < currSol.sequence.size() - 1 ; b++ ) 
		{
			int j = currSol.sequence[b];
			int l = currSol.sequence[b+1];
			double neighDecrement = - tsp.cost[h][i] - tsp.cost[j][l] + tsp.cost[h][j] + tsp.cost[i][l] ;
			if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) && !(neighDecrement < aspiration-0.01) ) 
			{
				continue;             
			} 
			if ( neighDecrement < bestDecrement ) 
			{
				bestDecrement = neighDecrement;
				move.from = a;
				move.to = b;
			}
		}
	}
	if ( bestDecrement >= tsp.infinite ) std::cout << "\n AARRGH!!! " << std::endl;
	return bestDecrement;
}

//Trova la minor differenza, non il vicino
double TSPSolver::findBestNeighbor3opt ( const TSP& tsp , const TSPSolution& currSol , TSPMove& move )
/* Determine the *move* yielding the best 3-opt neigbor solution
 */
{
	double bestDecrement = tsp.infinite;
	// intial and final position are fixed (initial/final node remains 0)
	for ( uint a = 1 ; a < currSol.sequence.size() - 3 ; a++ ) 
	{
		int h = currSol.sequence[a-1];
		int i = currSol.sequence[a];
		for ( uint b = a + 1 ; b < currSol.sequence.size() - 2 ; b++ ) 
		{
			int j = currSol.sequence[b];
			int k = currSol.sequence[b+1];
			//c parte da b+2 perché i 3 archi non devono essere consecutivi
			for ( uint c = b + 2 ; c < currSol.sequence.size() - 1 ; c++ ) 
			{
				int l = currSol.sequence[c];
				int m = currSol.sequence[c+1];
				//Nel 2opt ho solo un possibile scambio, ora ne ho 4 (ne avrei di più, ma se collego 2 lettere vicine ottengo sottocicli oppure posso avere dei 2opt che vanno scartati)
				double neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][l] + tsp.cost[h][j] + tsp.cost[m][k];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 1;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][m] + tsp.cost[l][j] + tsp.cost[h][k];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 2;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][k] + tsp.cost[h][l] + tsp.cost[m][j];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 3;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][l] + tsp.cost[h][k] + tsp.cost[m][j];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 4;
				}
			}
		}
	}
	return bestDecrement;
}

//Trova la minor differenza, non il vicino
double TSPSolver::findBestNeighbor3opt2 ( const TSP& tsp , const TSPSolution& currSol ,int currIter, TSPMove& move )
/* Determine the *move* yielding the best 3-opt neigbor solution
 */
{
	double bestDecrement = tsp.infinite;
	// intial and final position are fixed (initial/final node remains 0)
	for ( uint a = 1 ; a < currSol.sequence.size() - 3 ; a++ ) 
	{
		int h = currSol.sequence[a-1];
		int i = currSol.sequence[a];
		for ( uint b = a + 1 ; b < currSol.sequence.size() - 2 ; b++ ) 
		{
			int j = currSol.sequence[b];
			int k = currSol.sequence[b+1];
			//c parte da b+2 perché i 3 archi non devono essere consecutivi
			for ( uint c = b + 2 ; c < currSol.sequence.size() - 1 ; c++ ) 
			{
				int l = currSol.sequence[c];
				int m = currSol.sequence[c+1];
				//Nel 2opt ho solo un possibile scambio, ora ne ho 4 (ne avrei di più, ma se collego 2 lettere vicine ottengo sottocicli oppure posso avere dei 2opt che vanno scartati)
				if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) ) continue;
				double neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][l] + tsp.cost[h][j] + tsp.cost[m][k];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 1;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][m] + tsp.cost[l][j] + tsp.cost[h][k];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 2;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][k] + tsp.cost[h][l] + tsp.cost[m][j];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 3;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][l] + tsp.cost[h][k] + tsp.cost[m][j];
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 4;
				}
			}
		}
	}
	  
	return bestDecrement;
}

//Trova la minor differenza, non il vicino
double TSPSolver::findBestNeighbor3opt3 ( const TSP& tsp , const TSPSolution& currSol ,int currIter, double aspiration, TSPMove& move )
/* Determine the *move* yielding the best 3-opt neigbor solution
 */
{
	double bestDecrement = tsp.infinite;
	// intial and final position are fixed (initial/final node remains 0)
	for ( uint a = 1 ; a < currSol.sequence.size() - 3 ; a++ ) 
	{
		int h = currSol.sequence[a-1];
		int i = currSol.sequence[a];
		for ( uint b = a + 1 ; b < currSol.sequence.size() - 2 ; b++ ) 
		{
			int j = currSol.sequence[b];
			int k = currSol.sequence[b+1];
			//c parte da b+2 perché i 3 archi non devono essere consecutivi
			for ( uint c = b + 2 ; c < currSol.sequence.size() - 1 ; c++ ) 
			{
				int l = currSol.sequence[c];
				int m = currSol.sequence[c+1];
				//Nel 2opt ho solo un possibile scambio, ora ne ho 4 (ne avrei di più, ma se collego 2 lettere vicine ottengo sottocicli oppure posso avere dei 2opt che vanno scartati)
				double neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][l] + tsp.cost[h][j] + tsp.cost[m][k];
				if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) && !(neighDecrement < aspiration-0.01) ) 
				{
					continue;             
				} 
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 1;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][m] + tsp.cost[l][j] + tsp.cost[h][k];
				if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) && !(neighDecrement < aspiration-0.01) ) 
				{
					continue;             
				} 
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 2;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][k] + tsp.cost[h][l] + tsp.cost[m][j];
				if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) && !(neighDecrement < aspiration-0.01) ) 
				{
					continue;             
				} 
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 3;
				}
				neighDecrement = - tsp.cost[h][i] - tsp.cost[j][k] - tsp.cost[l][m] + tsp.cost[i][l] + tsp.cost[h][k] + tsp.cost[m][j];
				if ( (currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) && !(neighDecrement < aspiration-0.01) ) 
				{
					continue;             
				} 
				if ( neighDecrement < bestDecrement ) 
				{
					bestDecrement = neighDecrement;
					move.from = a;
					move.to = b;
					move.c = b+1;
					move.d=c;
					move.comb= 4;
				}
			}
		}
	}
	if ( bestDecrement >= tsp.infinite ) std::cout << "\n AARRGH!!! " << std::endl;
	return bestDecrement;
}

//3-opt. Scambia tre elementi non consecutivi della sequenza
TSPSolution& TSPSolver::swap3opt ( TSPSolution& tspSol , const TSPMove& move ) 
{
	//0 1 2 3 4 5 6 7 8 0
	//_ h i _ j k _ l m _
	//Le 4 combinazioni hanno degli swap particolari per ricreare il ciclo correttamente
	TSPSolution tmpSol(tspSol); 
	int comb=move.comb;
	if (comb==1)//i->j=j->i e k->l=l->k
	{
		cout<<endl<<"1) i->j=j->i e k->l=l->k"<<endl;
		tspSol.print();
		cout<<endl<<"i="<<move.from<<", j="<<move.to<<", k="<<move.c<<", l="<<move.d<<endl;
		for ( int i = move.from ; i <= move.to ; ++i ) 
		{
			tspSol.sequence[i] = tmpSol.sequence[move.to-(i-move.from)];
		}
		for ( int i = move.c ; i <= move.d ; ++i ) 
		{
			tspSol.sequence[i] = tmpSol.sequence[move.d-(i-move.c)];
		}
		tspSol.print();
	}
	else if (comb==2) //i->j=k->l e k->l=j->i
	{
		cout<<endl<<"2) i->j=k->l e k->l=j->i"<<endl;
		tspSol.print();
		cout<<endl<<"i="<<move.from<<", j="<<move.to<<", k="<<move.c<<", l="<<move.d<<endl;
		for(int i=0; i<=(move.d-move.c);++i)
		{
			tspSol.sequence[move.from+i]=tmpSol.sequence[move.c+i];
		}
		for(int i=0; i<=(move.to-move.from);++i)
		{
			tspSol.sequence[move.from+(move.d-move.c+1)+i]=tmpSol.sequence[move.to-i];
		}
		tspSol.print();
	}
	else if (comb==3) //i->j=l->k e k->l=i->j
	{
		cout<<endl<<"3) i->j=l->k e k->l=i->j"<<comb<<endl;
		tspSol.print();
		cout<<endl<<"i="<<move.from<<", j="<<move.to<<", k="<<move.c<<", l="<<move.d<<endl;
		for(int i=0; i<=(move.d-move.c);++i)
		{
			tspSol.sequence[move.from+i]=tmpSol.sequence[move.d-i];
		}
		for(int i=0; i<=(move.to-move.from);++i)
		{
			tspSol.sequence[move.from+(move.d-move.c+1)+i]=tmpSol.sequence[move.from+i];
		}
		tspSol.print();
	}
	else if (comb==4) //i->j=k->l e k->l=i->j
	{
		cout<<endl<<"4) i->j=k->l e k->l=i->j"<<endl;
		tspSol.print();
		cout<<endl<<"i="<<move.from<<", j="<<move.to<<", k="<<move.c<<", l="<<move.d<<endl;
		for(int i=0; i<=(move.d-move.c);++i)
		{
			tspSol.sequence[move.from+i]=tmpSol.sequence[move.c+i];
		}
		for(int i=0; i<=(move.to-move.from);++i)
		{
			tspSol.sequence[move.from+(move.d-move.c+1)+i]=tmpSol.sequence[move.from+i];
		}
		tspSol.print();
	}
	return tspSol;
}
