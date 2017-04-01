#include "info.h"
#include "parameter.h"

// 打印模拟时所需的各参数信息
void printPara(FILE* f, Parameter* para){

	if (! ifZeroRank())
      	return;

   	fprintf(f,
           "---Parameters information:---\n\n"
           "potentialName: %s\n"
           "xLatticeNum: %d\n"
           "yLatticeNum: %d\n"
           "zLatticeNum: %d\n"
           "xProcessNum: %d\n"
           "yProcessNum: %d\n"
           "zProcessNum: %d\n"
           "stepNums: %d\n"
           "printNums: %d\n"
           "stepTime: %g fs\n"
           "latticeConstant: %g\n"
           "initialTemperature: %g K\n"
           "initialDisplacement: %g Angstroms\n"
           "----------------\n",
           para->potentialName,
           para->xLat, 
           para->yLat,
           para->zLat,
           para->xProc,
           para->yProc,
           para->zProc,
           para->stepNums,
           para->printNums,
           para->stepTime,
           para->latConst,
           para->initTemper,
           para->initDisplace
   );
   fflush(f);

}