#include "atom.h"
#include "timer.h"
#include "parameter.h"
#include "cell.h"
#include "system.h"
#include "lattice.h"
#include "energy.h"

#include <stdlib.h>
#include <math.h>

#define KB (8.6173324e-5) //波尔兹曼常数

// 初始化原子信息结构体
void initAtoms(struct CellStr* cells, Atom** ato){

	*ato = (Atom*)malloc(sizeof(Atom));
    Atom* atoms = *ato;

   	int maxAtomNum = MAXPERCELL*cells->totalCellNum;
	
	atoms->myNum = 0;
   	atoms->totalNum = 0;

   	atoms->pos = (double3*) malloc(maxAtomNum*sizeof(double3));
   	atoms->momenta = (double3*) malloc(maxAtomNum*sizeof(double3));
   	atoms->force = (double3*) malloc(maxAtomNum*sizeof(double3));
   	atoms->pot = (double*)malloc(maxAtomNum*sizeof(double));
   	atoms->id = (int*)malloc(maxAtomNum*sizeof(int));

   	for (int i = 0; i < maxAtomNum; i++)
   	{
      	for(int j = 0; j< 3; j++){
      		atoms->pos[i][j] = 0.0;
      		atoms->momenta[i][j] = 0.0;
      		atoms->force[i][j] = 0.0;
      	}
      	atoms->pot[i] = 0.0;
      	atoms->id[i] = 0;
   	}
}

// 分配各原子到对应的细胞中
void distributeAtoms(struct SystemStr* sys, struct ParameterStr* para){
 
   	double latticeConst = sys->lattice->latticeConst;
   	int yLat = para->yLat;
   	int zLat = para->zLat;
   	double* myMin = sys->space->myMin;
   	double* myMax = sys->space->myMax;

   	double3 xyzpos;  // 原子坐标

   	double3 momenta; // 原子动量
   	momenta[0] = 0.0;
   	momenta[1] = 0.0;
   	momenta[2] = 0.0;
   
   	int n = 4;  // 每个晶胞4个原子
   	double3 displace[4] = { {0.25, 0.25, 0.25},
      	{0.25, 0.75, 0.75},
      	{0.75, 0.25, 0.75},
      	{0.75, 0.75, 0.25} };

   	// 分配原子
   	int begin[3];
   	int end[3];
   	for (int i = 0; i < 3; i++)
   	{
      	begin[i] = floor(myMin[i]/latticeConst);
      	end[i]   = ceil (myMax[i]/latticeConst);
   	}

   	for (int ix=begin[0]; ix<end[0]; ++ix)
      	for (int iy=begin[1]; iy<end[1]; ++iy)
         	for (int iz=begin[2]; iz<end[2]; ++iz)
            	for (int ib=0; ib<n; ++ib)
            	{
               		double xpos = (ix+displace[ib][0]) * latticeConst;
               		double ypos = (iy+displace[ib][1]) * latticeConst;
               		double zpos = (iz+displace[ib][2]) * latticeConst;
               		if (xpos < myMin[0] || xpos >= myMax[0]) continue;
               		if (ypos < myMin[1] || ypos >= myMax[1]) continue;
               		if (zpos < myMin[2] || zpos >= myMax[2]) continue;

               		// 计算原子的id
               		int id = ib+n*(iz+zLat*(iy+yLat*(ix)));

               		xyzpos[0] = xpos;
               		xyzpos[1] = ypos;
               		xyzpos[2] = zpos;

               		// 将此原子置于对应的细胞中,并初始化动量为0
               		assignAtom(id, xyzpos, sys, momenta);
            	}

   	// 利用mpi的reduce计算所有进程的总原子数量
   	//beginTimer(reduce);
   	//addIntParallel(&s->atoms->nLocal, &s->atoms->nGlobal, 1);
   	//endTimer(reduce);

   	//assert(s->atoms->nGlobal == nb*nx*ny*nz);
}

// 将指定原子分配到对应的细胞中
void assignAtom(int id, double3 xyzpos, struct SystemStr* sys, double3 momenta){
    
    // 根据原子坐标找到对应的细胞
    int cell = fineCellByCoord(sys->cells, sys->space, xyzpos);

    // 计算此原子为本空间第几个原子
    int n = cell*MAXPERCELL;
    n = n + sys->cells->atomNum[cell];
   
    // 若不在通信区域中，本空间总原子数加1
    if (cell < sys->cells->myCellNum)
        sys->atoms->myNum++;

    // 当前细胞中的原子数加1
    sys->cells->atomNum[cell]++;

    sys->atoms->id[n] = id;

    // 对原子的位置坐标、动量赋值
    for(int i =0; i<3 ;i++){
        sys->atoms->pos[n][i] = xyzpos[i];
        sys->atoms->momenta[n][i] = momenta[i];
    }
}

// 初始化体系的温度，即原子的速度
void initTemperature(struct SystemStr* sys, struct ParameterStr* para){

    // 指定温度
    double temper = para->initTemper;
    // 原子质量
    double atomM = sys->lat->atomM; 

    // 本空间所有原子总动量
    double3 myMomenta = {0.0,0.0,0.0};

    // 整个体系所有原子总动量
    double3 globalMomenta = {0.0,0.0,0.0};

    // 给定原子一个随机的速度及动量，并计算本空间的所有原子总动量
    for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
        for (int n=MAXPERCELL*nCell, count=0; count<sys->cells->atomNum[nCell]; count++, n++)
        {
            double sigma = sqrt(kB * temper/atomM);
            uint64_t seed = mkSeed(sys->atoms->id[n], 123);
            sys->atoms->momenta[n][0] = atomM * sigma * gasdev(&seed);
            sys->atoms->momenta[n][1] = atomM * sigma * gasdev(&seed);
            sys->atoms->momenta[n][2] = atomM * sigma * gasdev(&seed);

            myMomenta[0] += sys->atoms->momenta[n][0];
            myMomenta[1] += sys->atoms->momenta[n][1];
            myMomenta[2] += sys->atoms->momenta[n][2];
        }

    // 保证体系的总动量为0，在计算力之前需要调整为0

    // AllReduce, 得到整个体系的总动量
    MPI_Allreduce(myMomenta, globalMomenta, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // 每个原子需要中和的动量值
    double3 adjustMomenta;
    for(int i=0; i<3; i++)
        adjustMomenta[i] = -1 * globalMomenta[i]/sys->atoms->totalNum;

    // 调整各原子动量
    for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
        for (int n=MAXPERCELL*nCell, count=0; count<sys->cells->atomNum[nCell]; count++, n++)
            for(int i=0 ;i<3 ;i++)
                sys->atoms->momenta[n][i] += adjustMomenta[i];

    // 调整总动量为0后，需要调整体系的温度为指定温度
    computeTotalKinetic(sys);
    // 调整前的系统温度
    double t = (2*sys->energy->kineticEnergy)/(sys->atoms->totalNum*kB*3); 
    // 校正因子
    double factor = sqrt(temper/t);

    // 调整温度,乘以校正因子
    for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
        for (int n=MAXPERCELL*nCell, count=0; count<sys->cells->atomNum[nCell]; count++, n++)
            for(int i=0 ;i<3 ;i++)
                sys->atoms->momenta[n][i] *= factor; 

    // 计算调整后的总动能
    computeTotalKinetic(sys);
}

// 初始化原子的位移
void initDisplace(struct SystemStr* sys, struct ParameterStr* para){

}