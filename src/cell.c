#include "cell.h"

#include "space.h"
#include "potential.h"
#include <stdlib.h>

// 初始化细胞链表
void initCells(struct SpacialStr* space, struct PotentialStr* potential, struct CellStr* cells){

	cells = (Cell*)malloc(sizeof(Cell));

	// 保证细胞长度大于等于截断距离
	for (int i = 0; i < 3; i++)
   	{
      	cells->xyzCellNum[i] = space->myLength[i] / potential->cutoff; 
      	cells->cellLength[i] = space->myLength[i] / ((double) cells->xyzCellNum[i]);
   	}

   	// 实际细胞数为 x * y * z
   	cells->myCellNum = cells->xyzCellNum[0] * cells->xyzCellNum[1] * cells->xyzCellNum[2];
   
   	// 通信细胞数
   	cells->commCellNum = 2 * ((cells->xyzCellNum[0] + 2) *
                         (cells->xyzCellNum[1] + cells->xyzCellNum[2] + 2) +
                         (cells->xyzCellNum[1] * cells->xyzCellNum[2]));

   	// 总细胞数
   	cells->totalCellNum = cells->myCellNum + cells->commCellNum;
   
   	cells->atomNum = malloc(cells->totalCellNum*sizeof(int));

   	// 初始化各细胞中原子数为0
   	for (int i = 0; i < cells->totalCellNum; i++)
      	cells->atomNum[i] = 0;
}