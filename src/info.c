#include "info.h"
#include "mympi.h"
#include "atom.h"


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
           "initialTemperature: %g K\n"
           "----------------\n\n",
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
           para->initTemper
    );
    fflush(f);

}

// 打印势函数的相关信息
void printPotential(FILE* f, Potential* potential){

    if (! ifZeroRank())
        return;

    fprintf(f, "---Potential information:---\n\n");
    fprintf(f, "Potential type   : %s\n", potential->potentialType);
    fprintf(f, "Cutoff           : %g\n", potential->cutoff);
    fprintf(f, "sigma          : %g\n", potential->sigma);
    fprintf(f, "epsilon            : %g\n", potential->epsilon);
    //fprintf(f, "Beta            : %g\n", potential->Beta);
    fprintf(f, "----------------\n\n");
}

// 打印所模拟晶格的相关信息
void printLattice(FILE* f, Lattice* lattice){

    if (! ifZeroRank())
        return;

    fprintf(f, "---Lattice information:---\n\n");
    fprintf(f, "Lattice type    : %s\n", lattice->latticeType);
    fprintf(f, "Atom name       : %s\n", lattice->atomName);
    fprintf(f, "Atomic mass     : %g\n", lattice->atomM);
    fprintf(f, "Lattice Constant: %g\n", lattice->latticeConst);
    fprintf(f, "----------------\n\n");
}

// 跟踪模拟体系的总原子数
void printTotalAtom(FILE* f, Atom* atoms){
    if (! ifZeroRank())
        return;

    fprintf(f, "Total Atom    : %d\n", atoms->totalNum);
}

// 输出体系的温度
void printTemper(FILE*f, Energy* ener, int totalAtom){
    if (! ifZeroRank())
        return;

    double temper = (2*ener->kineticEnergy)/(totalAtom*kB*3);

    fprintf(f, "Temperature    : %g\n", temper);
}