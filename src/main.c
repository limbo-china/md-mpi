#include "timer.h"
#include "mympi.h"
#include "parameter.h"
#include "info.h"
#include "atom.h"
#include "potential.h"
#include "system.h"

#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

void updateMomenta(sys, para); 
void updatePosition(sys, para);

int main(int argc, char** argv){
	
	MPI_Init(&argc, &argv);
	initRank();

	//fprintf(stdout, "rankNums: %d\n", getRankNums());
	//fprintf(stdout, "myRank: %d\n", getMyRank());

	// 同步不起作用？
	//parallelBarrier("Begin to read parameters.");

	beginTimer(total);

	Parameter* para = readParameter();
	printPara(stdout,para);

	beginTimer(loop);
	//sleep(5);
	System* sys = initSystem(para);

	for(int i=0;i<10;i++){

    	updateMomenta(sys, para); 

    	updatePosition(sys, para);

    	adjustAtoms(sys);

    	computeForce(sys);

    	updateMomenta(sys, para); 

    	MPI_Allreduce(&sys->atoms->myNum, &sys->atoms->totalNum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    	printTotalAtom(stdout,sys->atoms);
    }
	endTimer(loop);
	

	endTimer(total);

	//fprintf(stdout, "total time: %g\n",getGlobalTime(total));
	//fprintf(stdout, "loop time: %g\n",getGlobalTime(loop));

	MPI_Finalize();
	return 0;
}

void updateMomenta(sys, para){

	double t = 0.5*para->stepTime;

	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0;i<3;i++)
         		sys->atoms->momenta[n][i] += t*sys->atoms->force[n][i];
}
void updatePosition(sys, para){

	double t = para->stepTime;
	double m = sys->lattice->atomM;

	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0;i<3;i++)
         		sys->atoms->pos[n][i] += t*sys->atoms->momenta[n][i]/m;
}