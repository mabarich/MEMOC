/**
 * @file TSPSolver.cpp
 * @brief TSP solver (neighborhood search)
 *
 */

#include "TSPSolver.h"
#include <iostream>


bool TSPSolver::solve ( const TSP& tsp , const TSPSolution& initSol , TSPSolution& bestSol )
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
		while ( ! stop ) 
		{
			if ( tsp.n < 20 ) currSol.print();
			std::cout << " (" << ++iter << ") value " << currValue << " (" << evaluate(currSol,tsp) << ")";
			//Passa ogni vicino e trova il miglior risultato
			double bestNeighValue = currValue + findBestNeighbor(tsp,currSol,move);
			std::cout << " move: " << move.from << " , " << move.to << std::endl;
			//Criterio di stop: se ho un risultato migliorante, vado la, altrimenti mi fermo
			if ( bestNeighValue < currValue )
			{
				bestValue = currValue = bestNeighValue;
				//Passo al vicino
				currSol = swap(currSol,move);
				stop = false;
			} 
			else 
			{
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

bool TSPSolver::solveTL ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol )
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
}

bool TSPSolver::solveAS ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol )   
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
