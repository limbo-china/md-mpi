//timer.h 
//作为程序各计算部分的计时器，用于性能分析

#ifndef TIMER_H_
#define TIMER_H_
	
//计时器数组指针
enum TimerPtr{
	total,
	loop,
	computeForce,
	atomExchange,
	reduce,
	nums
};

//计时器结构体,单位均为us
typedef struct TimerStr{
	uint64_t begin; //开始时间
	uint64_t global; //所有调用的总时间消耗
	uint64_t delta; //最近一次调用的时间消耗
}Timer;

void beginTimer(const enum TimerPtr ptr);
void endTimer(const enum TimerPtr ptr);

#endif