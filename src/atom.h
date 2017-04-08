// atom.h
// 原子结构体，存储原子的作用力，势能，速度等信息

#ifndef ATOM_H_
#define ATOM_H_

// 每个细胞中的原子数最大值
#define MAXPERCELL 64

#include "mytype.h"

struct CellStr;
struct SystemStr;
struct ParameterStr;

typedef struct AtomStr{

	double3*  pos;     // 原子坐标
   	double3*  momenta;     // 原子动量
   	double3*  force;     // 原子受到的作用力 
   	double*  pot;     // 原子势能

	int myNum; // 本进程空间中的总原子数
	int totalNum; // 整个体系的总原子数

	int* id;      // 各原子id

}Atom;

// 初始化原子信息
void initAtoms(struct CellStr* cells, Atom** ato);

// 分配各原子到对应的细胞中
void distributeAtoms(struct SystemStr* sys, struct ParameterStr* para);

// 将指定原子根据其坐标，分配到对应的细胞中
void assignAtom(int id, double3 xyzpos, struct SystemStr* sys, double3 momenta);

#endif