#include "timer.h"
#include "mympi.h"
#include "parameter.h"
#include "info.h"

#include <stdio.h>
#include <time.h>

int main(int* argc, char** argv){
	
	MPI_Init(&argc, &argv);
	initRank();

	fprintf(stdout, "rankNums: %d\n", getRankNums());
	fprintf(stdout, "myRank: %d\n", getMyRank());

	beginTimer(total);
	
	Parameter* para = parseParameter();

	beginTimer(loop);
	sleep(5);
	endTimer(loop);

	printPara(stdout,para);

	endTimer(total);

	fprintf(stdout, "total time: %g\n",getGlobalTime(total));
	fprintf(stdout, "loop time: %g\n",getGlobalTime(loop));
	return 0;
}