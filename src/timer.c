#include "timer.h"

#include <stdio.h>
#include <time.h>
#include <string.h>

static Timer timers[nums]; //计时器数组

static uint64_t getUsTime(); //获取当前时间(us)
static double usecToSec(uint64_t t); //us转化为s

static uint64_t getUsTime(){
	struct timeval t;
   	gettimeofday(&t, (struct timezone *)NULL);
   	return ((uint64_t)1000000)*(uint64_t)t.tv_sec + (uint64_t)t.tv_usec; 
}

static double usecToSec(uint64_t t){
	return t*1.0e-6;
}

void beginTimer(const enum TimerPtr ptr){
	timers[ptr].begin = getUsTime();
}

void endTimer(const enum TimerPtr ptr){
	timers[ptr].delta = getUsTime() - timers[ptr].begin;
	timers[ptr].global = timers[ptr].global + timers[ptr].delta;
}