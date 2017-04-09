#include "datacomm.h"
#include "cell.h"

#include <stdlib.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

// 初始化结构体
void initComm(DataComm** comm, struct SpacialStr* space, struct CellStr* cells){

	*comm = (DataComm*)malloc(sizeof(DataComm));
    DataComm* datacomm = *comm;

    int* myPos = space->position;
    int* globalProcNum = space->globalProcNum;
    int* xyzCellNum = cells->xyzCellNum;

    // 计算邻居进程号
    datacomm->neighborProc[X_NEG] = (myPos[0] -1 + globalProcNum[0]) % globalProcNum[0]
    	+ globalProcNum[0] *(myPos[1] + globalProcNum[1]*myPos[2]);
    datacomm->neighborProc[X_POS] = (myPos[0] +1 + globalProcNum[0]) % globalProcNum[0]
    	+ globalProcNum[0] *(myPos[1] + globalProcNum[1]*myPos[2]);

    datacomm->neighborProc[Y_NEG] = myPos[0] + globalProcNum[0] *
    	( (myPos[1] -1 + globalProcNum[1]) % globalProcNum[1] + globalProcNum[1]*myPos[2]);
    datacomm->neighborProc[Y_POS] = myPos[0] + globalProcNum[0] *
    	( (myPos[1] +1 + globalProcNum[1]) % globalProcNum[1] + globalProcNum[1]*myPos[2]);

    datacomm->neighborProc[Z_NEG] = myPos[0] + globalProcNum[0] *
    	( myPos[1] + globalProcNum[1]*((myPos[2] -1 + globalProcNum[2]) % globalProcNum[2]));
    datacomm->neighborProc[Z_POS] = myPos[0] + globalProcNum[0] *
    	( myPos[1] + globalProcNum[1]*((myPos[2] +1 + globalProcNum[2]) % globalProcNum[2]));

    // 各方向需要通信的细胞数的最大值
    int maxComm = MAX((xyzCellNum[0]+2)*(xyzCellNum[1]+2),
    	MAX((xyzCellNum[1]+2)*(xyzCellNum[2]+2),
    		(xyzCellNum[0]+2)*(xyzCellNum[2]+2)));
    datacomm->bufSize = 2*maxComm*MAXPERCELL*sizeof(AtomData);

    datacomm->commCellNum[X_NEG] = 2*(xyzCellNum[1]+2)*(xyzCellNum[2]+2);
   	datacomm->commCellNum[X_POS] = 2*(xyzCellNum[1]+2)*(xyzCellNum[2]+2);
   	datacomm->commCellNum[Y_NEG] = 2*(xyzCellNum[0]+2)*(xyzCellNum[2]+2);
   	datacomm->commCellNum[Y_POS]  = 2*(xyzCellNum[0]+2)*(xyzCellNum[2]+2);
   	datacomm->commCellNum[Z_NEG]  = 2*(xyzCellNum[0]+2)*(xyzCellNum[1]+2);
   	datacomm->commCellNum[Z_POS]  = 2*(xyzCellNum[0]+2)*(xyzCellNum[1]+2);

   	for (int dimen=0; dimen<6; dimen++)
      datacomm->commCells[dimen] = findCommCells(cells, dimen, datacomm->commCellNum[dimen]);
}

// 找出指定维度上所有通信部分的细胞
int* findCommCells(struct CellStr* cells, enum Neighbor dimen, int num){
	
	int* commcells = malloc(num*sizeof(int));
   	int xBegin = -1;
   	int xEnd   = cells->xyzCellNum[0]+1;
   	int yBegin = -1;
   	int yEnd   = boxes->xyzCellNum[1]+1;
   	int zBegin = -1;
   	int zEnd   = boxes->xyzCellNum[2]+1;

   	if (dimen == X_NEG) xEnd = xBegin+2;
   	if (dimen == X_POS) xBegin = xEnd-2;
   	if (dimen == Y_NEG) yEnd = yBegin+2;
   	if (dimen == Y_POS) yBegin = yEnd-2;
   	if (dimen == Z_NEG) zEnd = zBegin+2;
   	if (dimen == Z_POS) zBegin = zEnd-2;

   	int n = 0;
   	for (int ix=xBegin; ix<xEnd; ix++)
      	for (int iy=yBegin; iy<yEnd; iy++)
         	for (int iz=zBegin; iz<zEnd; iz++)
            	commcells[n++] = findCellByXYZ(cells, ix, iy, iz);
   	//assert
   	return commcells;
}
