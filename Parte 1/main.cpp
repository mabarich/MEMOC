/**
 * @file antenne.cpp
 * @brief 
 */

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <fstream>
#include <math.h> 
#include "cpxmacro.h"

using namespace std;

// error status and messagge buffer
int status;
char errmsg[BUF_SIZE];
int fun;
//Numero nodi
//const int N = 50; 
//Numero di cluster  
int CL=1;   
//Etichette nodi
//double x[N];
//double y[N];
//double dist[(N*N)];

int N; //number of nodes
std::vector<double> x;
std::vector<double> y;
std::vector<double> dist;

//Nomi variabili			
const int NAME_SIZE = 512;
char name[NAME_SIZE];

void distances()
{
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			dist[(i*N)+j]=sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
		}
	}
}

void writePoints()
{
	ofstream output_file;
	output_file.open("output.txt",ios::out);
	for(int i=0; i<N; i++)
	{
		output_file<<x[i]<< " "<<y[i]<<endl;
	}
}
	
void createPoints()
{
	srand ( time(NULL) );
	int min=0;
	int max=10;
	for(int i=0; i<N; i++)
	{
		x[i]=rand()%(max-min + 1) + min;
	} 
	for(int i=0; i<N; i++)
	{
		y[i]=rand()%(max-min + 1) + min;
	} 
	distances();
}	

void createClPoints()
{
	int cluster=0;
	int xcl[CL];
	int ycl[CL];
	//Creo un punto come centro di ogni cluster (11,11), (22,22), ecc...
	for(int i=0; i<CL; i++)
	{
		xcl[i]=11*(i+1);
		ycl[i]=11*(i+1);
	}
	int pos=0;
	int resto=N%CL;
	int max=4;
	int min=0;
	while(cluster<CL)
	{
		srand ( time(NULL) );
		int cont=0;
		while (cont<(N/CL))
		{
			double alpha = 2.0d*M_PI*((double) rand() / RAND_MAX);
			double radius=rand()%(max-min + 1) + min;
			x[pos] = radius*cos(alpha)+xcl[cluster];
			y[pos] = radius*sin(alpha)+ycl[cluster];
			pos++;
			cont++;
		}
		cluster++;
	}
	cluster--;
	for (int i=0; i<resto; i++)
	{
		double alpha = 2.0d*M_PI*((double) rand() / RAND_MAX);
		double radius=rand()%(max-min + 1) + min;
		x[pos] = radius*cos(alpha)+xcl[cluster];
		y[pos] = radius*sin(alpha)+ycl[cluster];
		pos++;
	}
	distances();
}

void create4Points()
{
	for(int i=0; i<4; i++)
	{
		x[i]=i;
	} 
	y[0]=1;
	y[1]=8;
	y[2]=5;
	y[3]=11;
	distances();
}

void createCPoints()
{
	srand ( time(NULL) );
	int xcenter=3;
	int ycenter=3;
	int radius=100;
	for(int i=0; i<N; i++)
	{
		double alpha = 2.0d*M_PI*((double) rand() / RAND_MAX);
		x[i] = radius*cos(alpha)+xcenter;
		y[i] = radius*sin(alpha)+ycenter;
	}
	distances();
}

void readFile()
{
	std::ifstream in("points.txt");
	in>>N;
	dist.reserve(N*N);
	x.reserve(N);
	y.reserve(N);
	for (int i = 0; i < N; i++) 
	{
		double xx;
		in >> xx;
		double yy;
		in >> yy;
		x[i]=xx;
		y[i]=yy;
	}
	in.close();
}
	
void setupLP(CEnv env, Prob lp, int & numVars )
{
	//Creo i punti e calcolo le distanze
	if (fun<=0)
		create4Points();
	else if (fun==1)
		createPoints();
	else if (fun==2)
		createCPoints();
	else if (fun==3)
		createClPoints();
	else if (fun>=4)
		readFile();
	writePoints();
	
	//Variabili xij dei flussi (in ordine)
	for (int n1 = 0; n1 < N; n1++)
	{
		for (int n2 = 0; n2 < N; n2++)
		{
			char xtype = 'C';
			double obj = 0.0;
			double lb = 0.0;
			double ub = N;
			snprintf(name, NAME_SIZE, "x_%d_%d", n1, n2);
			char* xname = (char*)(&name[0]);
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1, &obj, &lb, &ub, &xtype, &xname );
		}
	}
	//Variabili yij degli archi (in ordine) e funzione obiettivo
	for (int a1 = 0; a1 < N; a1++)
	{
		for (int a2 = 0; a2 < N; a2++)
		{
			char ytype = 'B';	
			double obj = 0.0;
			if(a1!=a2)
				obj = dist[a1*N+a2];
			double lb = 0.0;
			double ub = 1.0;
			snprintf(name, NAME_SIZE, "y_%d_%d", a1, a2);
			char* yname = (char*)(&name[0]);
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1, &obj, &lb, &ub, &ytype, &yname );
		}
	}
	numVars = CPXgetnumcols(env, lp);
	//Somm x0j=N
	//So che le x0j sono le prime N variabili
	{
		std::vector<int> idx;
		std::vector<double> coef;
		for (int i = 0; i < N; i++)
		{	
			idx.push_back(i);  
			coef.push_back(1);
		}
		char sense = 'E';
		double rhs = N;
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
	}
	//Somm xik-Som xkj=1
	for (int k = 1; k < N; k++)
	{
		std::vector<int> idx;
		std::vector<double> coef;
		//Sto esaminando i nodi uscenti da k
		for (int i = 0; i < N; i++)
		{	
			if(k!=i)
			{
				idx.push_back((k*N)+i);  
				coef.push_back(-1);
			}		
		}
		//Sto esaminando i nodi entranti in k
		for (int i = 0; i < N; i++)
		{	
			if(k!=i)
			{
				idx.push_back(k+(i*N));  
				coef.push_back(1);
			}
		}
		char sense = 'E';
		double rhs = 1;
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
	}
	//Somm yij=1 (per ogni j)
	for (int j = 0; j < N; j++)
	{
		std::vector<int> idx;
		std::vector<double> coef;
		for (int i = 0; i < N; i++)
		{	
			if(i!=j)
			{
				idx.push_back((N*N)+(j+(i*N))); 
				coef.push_back(1);
			}
		}
		char sense = 'E';
		double rhs = 1;
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
	}
	//Somm yij=1 (per ogni i)
	//So che le x sono ordinate, quindi non mi serve il vector ord ma basta pos
	for (int i = 0; i < N; i++)
	{
		std::vector<int> idx;
		std::vector<double> coef;
		for (int j = 0; j < N; j++)
		{	
			if(i!=j)
			{
				idx.push_back((N*N)+((i*N)+j)); 
				coef.push_back(1);
			}
		}
		char sense = 'E';
		double rhs = 1;
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
	}
	//xij<= N*yij
	//So che le x sono ordinate, quindi non mi serve il vector ord ma basta pos
	for (int n1 = 0; n1 < N; n1++)
	{
		for (int n2 = 0; n2 < N; n2++)
		{
			if(n1!=n2)
			{
				std::vector<int> idx;
				std::vector<double> coef;
				idx.push_back((n1*N)+n2);  
				coef.push_back(1);
				idx.push_back((N*N)+((n1*N)+n2));  
				coef.push_back((-1)*N);
				char sense = 'L';
				double rhs = 0.0;
				int matbeg = 0;
				CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
			}
		}
	}

	CPXchgobjsen(env, lp, CPX_MIN);
	// print (debug)
	CHECKED_CPX_CALL( CPXwriteprob, env, lp, "Esercitazione.lp", 0 );
}

int main (int argc, char const *argv[])
{
	try
	{
		if (argc < 2) throw std::runtime_error("usage: ./main #points #function [#clusters]");
		N=atoi(argv[1]);
		fun = atoi(argv[2]);
		if (fun==3)
		{
			CL=atoi(argv[3]);
		}
		if(fun!=4)
		{
			dist.reserve(N*N);
			x.reserve(N);
			y.reserve(N);
		}
		// init
		DECL_ENV( env );
		DECL_PROB( env, lp );
		// setup LP
		int numVars;
		setupLP(env, lp, numVars);
		// optimize
		CHECKED_CPX_CALL( CPXmipopt, env, lp );
		// print
		double objval;
		CHECKED_CPX_CALL( CPXgetobjval, env, lp, &objval );
		std::cout << "Objval: " << objval << std::endl;
		int n = CPXgetnumcols(env, lp);
		cout << n << " " << numVars << endl;
		if (n != numVars) { throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + "different number of variables"); }
	  std::vector<double> varVals;
	  varVals.resize(n);
  	CHECKED_CPX_CALL( CPXgetx, env, lp, &varVals[0], 0, n - 1 );
  	for ( int i = 0 ; i < n ; ++i ) {
  	  std::cout << "var in position " << i << " : " << varVals[i] << std::endl;
  	}
		CHECKED_CPX_CALL( CPXsolwrite, env, lp, "Esercitazione.sol" );
		// free
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}
