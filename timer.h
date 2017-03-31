#ifndef TIMER_H_
#define TIMER_H_
	
enum TimerPtr{
	global,
	loop,
	computeForce,
	atomExchange,
	reduce,
	nums
};

void startTimer(const enum TimerPtr ptr);
void stopTimer(const enum TimerPtr ptr);

#endif