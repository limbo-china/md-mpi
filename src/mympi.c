#include "mympi.h"

// 当前进程编号
static int myRank = 0;
// 并行总进程数
static int rankNums = 1;

// 获取并行的总进程数
int getRankNums()
{
   	return rankNums;
}

// 获取当前进程编号
int getMyRank()   
{
   	return myRank;
}

// 判断当前进程是否为0进程
// 用于控制台输出，只需要0进程输出即可
int ifZeroRank()
{
   	if (myRank == 0) 
   		return 1;
   	return 0;
}

//获取进程数和编号
void initRank()
{
   	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   	MPI_Comm_size(MPI_COMM_WORLD, &rankNums);
}