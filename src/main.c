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

void updateMomenta(System* sys, Parameter* para); 
void updatePosition(System* sys, Parameter* para);

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

	int calls=0;
	for(int i=0;i<para->stepNums;i++){
		calls++;
    	updateMomenta(sys, para); 

    	updatePosition(sys, para);

    	beginTimer(adjustatom);
    	adjustAtoms(sys);
    	endTimer(adjustatom);

    	//beginTimer(force);
    	computeForce(sys);
		//endTimer(force);

    	updateMomenta(sys, para); 
    	if(i%para->printNums == 0){

    	//MPI_Allreduce(&sys->atoms->myNum, &sys->atoms->totalNum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    	//printTotalAtom(stdout,sys->atoms);

    	computeTotalKinetic(sys);
    	printTemper(stdout,sys->energy,sys->atoms->totalNum);
    	}
    }
    printf("calls:%d",calls);
	endTimer(loop);
	

	endTimer(total);

	//fprintf(stdout, "total time: %g\n",getGlobalTime(total));
	fprintf(stdout, "loop time: %g\n",getGlobalTime(loop));
	fprintf(stdout, "adjust time: %g\n",getGlobalTime(adjustatom));
	fprintf(stdout, "comm time: %g\n",getGlobalTime(communication));
	fprintf(stdout, "force time: %g\n",getGlobalTime(force));



	MPI_Finalize();
	return 0;
}

void updateMomenta(System* sys, Parameter* para){

	double t = 0.5*para->stepTime;

	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0;i<3;i++)
         		sys->atoms->momenta[n][i] += t*sys->atoms->force[n][i];
}
void updatePosition(System* sys, Parameter* para){

	double t = para->stepTime;
	double m = sys->lattice->atomM;

	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0;i<3;i++)
         		sys->atoms->pos[n][i] += t*sys->atoms->momenta[n][i]/m;
}