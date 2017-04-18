#include "timer.h"
#include "mympi.h"
#include "parameter.h"
#include "info.h"
#include "system.h"

#include <stdio.h>
#include <unistd.h>

int main(int argc, char** argv){
	
	MPI_Init(&argc, &argv);
	initRank();

	fprintf(stdout, "rankNums: %d\n", getRankNums());
	fprintf(stdout, "myRank: %d\n", getMyRank());

	// 同步不起作用？
	//parallelBarrier("Begin to read parameters.");

	beginTimer(total);

	Parameter* para = readParameter();
	printPara(stdout,para);

	beginTimer(loop);
	//sleep(5);
	System* sys = initSystem(para);
	endTimer(loop);
	

	endTimer(total);

	//fprintf(stdout, "total time: %g\n",getGlobalTime(total));
	//fprintf(stdout, "loop time: %g\n",getGlobalTime(loop));

	MPI_Finalize();
	return 0;
}