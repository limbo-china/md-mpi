// mpi.h
// 程序并行的一些基本函数

#ifndef MPI_H_
#define MPI_H_

// 获取并行的总进程数
int getRankNums();

// 获取当前进程编号
int getMyRank();

// 判断当前进程是否为0进程
int ifZeroRank();

//获取进程数和编号
void initRank();
#endif