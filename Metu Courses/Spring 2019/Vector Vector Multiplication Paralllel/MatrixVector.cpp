#include<iostream>
#include <fstream>
#include <string>
# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <mpi.h>
//#include "SCmathlib.h"
//# include "blaswrap.h"
//# include "f2c.h"
#include <lapacke.h>
#include <cblas.h>


using namespace std;

int main(int argc, char *argv[]){
	int rank, buf, totalnodes;
	int startval,endval,accum;
	int i,j;
	// the size of array 
    int n;
	std::clock_t start;
    double duration;
	double sum;
	
	//cout<<"Enter the size of vector : "<<endl;
	//Taking input in array  
    //cin>>n; 
	n=4;
	// the vector to be computed
	double r[n];
	// the vector with n size
     for(int j=0; j<n; j++)
     {
        r[j]=1.*(j+1); 
     }
	
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // get totalnodes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get rank
		
	printf("Proc %d of %d > Does anyone have a toothpick?\n",rank, totalnodes);
	// let's output the total number of nodes
	cout << "The total number of nodes is: " << totalnodes << "\n" <<endl;
	cout << "The rank of nodes is: " << rank << "\n" <<endl;
	
	
	//sum = 0; // zero sum for accumulation
	//startval = 1000*rank/totalnodes+1;
	//endval = 1000*(rank+1)/totalnodes;

	double c =  cblas_ddot(n/totalnodes, r, 1, r, 1);
	MPI_Reduce(&c, &sum, 4, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Finalize();
	
	cout<<"the result of dot product "<< sum <<"\n"<<endl; 
	//for(int i=startval;i<=endval;i=i+1)
	//j=i;
	//sum = sum + i;
	//if(rank!=0)
		//cout << "My rank : " << rank << "\n" <<endl;
	//	MPI_Send(&sum,10,MPI_INT,0,1,MPI_COMM_WORLD);
	//else
	//	cout << "My rank : " << rank << "\n" <<endl;
	//	for(int j=1;j<totalnodes;j=j+1){			
	//		MPI_Recv(&accum,10,MPI_INT,j,1,MPI_COMM_WORLD, &status);
			//sum = sum + accum;
	//	}
	
	//MPI_Reduce(&c, &sum, 4, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
	//if(rank == 0)
	//cout << "The sum from 1 to 1000 is: " << sum << endl;
	
	
	
	
	
    
	// let's start the times
    //start = std::clock();
	// dot product of a vector 	
	//for(int i=startval;i<=endval;i=i+1){
	//	c =  cblas_ddot(n, r, 1, r, 1);
	//}
	//cout<<"the result of dot product "<<c <<"\n"<<endl; 
	// the  time elapsed during the dot product operation
	//duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout<<"seconds elapsed during the operation is: "<< duration <<'\n';
	
	// let's output the time
    //ofstream myfile;
  	//myfile.open ("SecondsPerCores.txt");
  	//myfile << "# of Cores" << duration <<' seconds \n';
  	//myfile.close();

	
return 0;
}