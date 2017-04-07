// atom.h
// 原子结构体，存储原子的作用力，势能，速度等信息

#ifndef ATOM_H_
#define ATOM_H_

#include "mytype.h"
#include "cell.h"
#include "system.h"
#include <stdlib.h>

// 每个细胞中的原子数最大值
#define MAXPERCELL 64

typedef struct AtomStr{

	double3*  pos;     // 原子坐标
   	double3*  momenta;     // 原子动量
   	double3*  force;     // 原子受到的作用力 
   	double3*  pot;     // 原子势能

	int myNum; // 本进程空间中的总原子数
	int totalNum; // 整个体系的总原子数

	int* id;      // 各原子id

}Atom;

// 初始化原子信息
void initAtoms(Cell* cells, Atom* atoms);

// 分配各原子到对应的细胞中
void distributeAtoms(System* sys, Parameter* para);

#endif