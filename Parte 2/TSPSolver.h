/**
 * @file TSPSolver.h
 * @brief TSP solver (neighborhood search)
 *
 */

#ifndef TSPSOLVER_H
#define TSPSOLVER_H

#include <vector>

#include "TSPSolution.h"

/**
 * Class representing substring reversal move
 */
//Mossa. Uso mosse e non nodi e basta così, anche se ho una mossa nella tabu list, posso passare per quei nodi per prendere altri vicini
typedef struct move 
{
  int			from;
  int			to;
} TSPMove;

/**
 * Class that solves a TSP problem by neighbourdood search and 2-opt moves
 */
class TSPSolver
{
public:
	/** Constructor */
	TSPSolver ( ) { }
	/**
	* evaluate a solution
	* @param solution solution to be evaluated
	* @param TSP TSP data
	* @return the value of the sequence
	*/
	//Somma i pesi degli archi
	double evaluate ( const TSPSolution& sol , const TSP& tsp ) const 
	{
		double total = 0.0;
		uint cont=sol.sequence.size();
		for ( uint i = 0 ; i < cont - 1 ; ++i ) 
		{
			int from = sol.sequence[i]  ;
			int to   = sol.sequence[i+1];
			total += tsp.cost[from][to];
		}
		return total;
	}
	/**
	* initialize a solution as a random sequence by random swaps
	* @param sol solution to be initialized
	* @return true if everything OK, false otherwise
	*/
	//Soluzione iniziale casuale
	bool initRnd ( TSPSolution& sol ) 
	{
		srand(time(NULL));
		for ( uint i = 1 ; i < sol.sequence.size() ; ++i ) 
		{
			// intial and final position are fixed (initial/final node remains 0)
			int idx1 = rand() % (sol.sequence.size()-2) + 1;
			int idx2 = rand() % (sol.sequence.size()-2) + 1;
			int tmp = sol.sequence[idx1];
			sol.sequence[idx1] = sol.sequence[idx2];
			sol.sequence[idx2] = tmp;
		}
		std::cout << "### "; sol.print(); std::cout << " ###" << std::endl;
		return true;
	}
	//Prendo sempre il migliore dei nodi rimanenti
	bool initBst ( TSPSolution& sol , const TSP& tsp ) 
	{
		vector<int> rim;
		rim.reserve(tsp.n-1);
		for(int i=1; i<tsp.n;i++)
		{
			rim.push_back(i);		
		}
		uint i=1;
		while(i<=tsp.n)
		{
			int bestnode=-1;
			double bestcost=tsp.infinite;
			for(uint j=0; j<rim.size(); j++)
			{
				double costo=tsp.cost[sol.sequence[i-1]][j];
				if(costo<=bestcost)
				{
					bestcost=costo;
					bestnode=j;	
				}			
			}		
			sol.sequence[i]=rim[bestnode];
			rim.erase(rim.begin()+bestnode);
			i++;
		}
		std::cout << "### "; sol.print(); std::cout << " ###" << std::endl;
		return true;
	}
	/**
	* Main method: search for a good tour by local search
	* @param TSP TSP data
	* @param initSol initial solution
	* @param bestSol best found solution (output)
	* @return true id everything OK, false otherwise
	*/
	bool solve ( const TSP& tsp , const TSPSolution& initSol , TSPSolution& bestSol );
	bool solveTL ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol );
	bool solveAS ( const TSP& tsp , const TSPSolution& initSol , int tabulength , int maxIter , TSPSolution& bestSol );
protected:
	double		findBestNeighbor ( const TSP& tsp , const TSPSolution& currSol , TSPMove& move );
	double		findBestNeighbor2 ( const TSP& tsp , const TSPSolution& currSol , int currIter , TSPMove& move );
	double		findBestNeighbor3 ( const TSP& tsp , const TSPSolution& currSol , int currIter , double aspiration , TSPMove& move );
	TSPSolution&  swap 			 ( TSPSolution& tspSol , const TSPMove& move );
	bool again(TSPMove& move, std::vector<TSPMove>& moves);
	
	///Tabu search (tabu list stores when a node is involved)
	int               tabuLength;
	std::vector<int>  tabuList;
	bool  initTabuList ( int n ) 
	{
		for ( int i = 0 ; i < n ; ++i ) 
		{
			tabuList.push_back(-tabuLength-1);
		}
		return true;
	}
	//Memoria per la diversificazione
	static const uint memLength=10;
};

#endif /* TSPSOLVER_H */
