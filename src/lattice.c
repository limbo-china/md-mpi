#include "lattice.h"

#include <stdlib.h>

// 初始化晶格结构体
void initLatticeInfo(Lattice* lattice){

	lattice = (Lattice*)malloc(sizeof(Lattice));

	strcpy(lattice->latticeType, "FCC");
	strcpy(lattice->atomName, "Cu");
	lattice->atomM = 63.55;
	lattice->latticeConst = 3.615;
}