#include "system.h"

#include <stdlib.h>
#include <stdio.h>

//初始化模拟体系
System* initSystem(Parameter* para){

	System* sys = (System*)malloc(sizeof(System));
	//sys->energy = NULL;
	sys->potential = NULL;
	sys->lattice = NULL;
    sys->space = NULL;
    sys->cells = NULL;
    sys->atoms = NULL;
    //sys->atomExchange = NULL;

   	initPotInfo(&sys->potential);  // 传值问题！！！！！！！！！
    printPotential(stdout, sys->potential);
   	initLatticeInfo(&sys->lattice);
    printLattice(stdout, sys->lattice);
    initSpace(para, sys->lattice, &sys->space);
    initCells(sys->space, sys->potential, &sys->cells);
    initAtoms(sys->cells, &sys->atoms);

    distributeAtoms(sys, para);

    return sys;
}