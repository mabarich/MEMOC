/**
 * @file TSP.h
 * @brief TSP data
 *
 */

#ifndef TSP_H
#define TSP_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h> 

using namespace std;

/**
 * Class that describes a TSP instance (a cost matrix, nodes are identified by integer 0 ... n-1)
 */
class TSP
{
public:
	/*static const int n=10; //number of nodes
	double x[n];
	double y[n];*/
	int n; //number of nodes
	std::vector<double> x;
	std::vector<double> y;
	std::vector< std::vector<double> > cost;
	void read(const int num, const int fun, const int cluster)
	{
		n=num;
		if(fun!=4)
		{
			cost.reserve(num);
			x.reserve(num);
			y.reserve(num);
		}
		if (fun<=0)
			create4Points();
		else if (fun==1)
			createPoints();
		else if (fun==2)
			createCPoints();
		else if (fun==3)
			createClPoints(cluster);
		else if (fun>=4)
			readFile();
		writePoints();
		distances();
	}

//Mette le distanze per ogni coppia di punti nel vector  
void distances()
{
	cost.resize(n);
	for(int i=0; i<n; i++)
	{
		cost[i].reserve(n);
		for(int j=0; j<n; j++)
		{
			cost[i].push_back(sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i])));
		}
	}
}

void writePoints()
{
	ofstream output_file;
	output_file.open("output.txt",ios::out);
	for(int i=0; i<n; i++)
	{
		output_file<<x[i]<< " "<<y[i]<<endl;
	}
}
	
//Punti random
void createPoints()
{
	srand ( time(NULL) );
	int min=0;
	int max=200;
	for(int i=0; i<n; i++)
	{
		x[i]=rand()%(max-min + 1) + min;
	} 
	for(int i=0; i<n; i++)
	{
		y[i]=rand()%(max-min + 1) + min;
	} 
}	

//Punti random in CL clusters (coordinate polari)
void createClPoints(const int CL)
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
	int resto=n%CL;
	int max=4;
	int min=0;
	while(cluster<CL)
	{
		srand ( time(NULL) );
		int cont=0;
		while (cont<(n/CL))
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
}

//Punti random lungo una circonferenza (coordinate polari)
void createCPoints()
{
	srand ( time(NULL) );
	int xcenter=3;
	int ycenter=3;
	int radius=3;
	for(int i=0; i<n; i++)
	{
		double alpha = 2.0d*M_PI*((double) rand() / RAND_MAX);
		x[i] = radius*cos(alpha)+xcenter;
		y[i] = radius*sin(alpha)+ycenter;
	}
}

void readFile()
{
	std::ifstream in("points.txt");
	in>>n;
	x.reserve(n);
	y.reserve(n);
	for (int i = 0; i < n; i++) 
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
  static const double infinite=1e10;
};

#endif /* TSP_H */

