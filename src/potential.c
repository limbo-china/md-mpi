#include "potential.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// 初始化势函数结构体
void initPotInfo(Potential** pot){


	*pot = (Potential*)malloc(sizeof(Potential));
	Potential* potential= *pot;
	
	strcpy(potential->potentialType,"Lennard-Jones");
	potential->sigma = 2.315;	                  
   	potential->epsilon = 0.167;
   	potential->cutoff = 2.5*potential->sigma;

   	//potential->computeforce = computeForce;
   	//potential->free = potentialFree;
}

// 释放结构体空间
void potentialFree(Potential* potential){
	if(potential)
		free(potential);
}