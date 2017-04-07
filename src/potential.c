#include "potential.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// 初始化势函数结构体
void initPotInfo(Potential* potential){

	fprintf(stdout, "Potential size : %d\n",sizeof(Potential) );
	potential = (Potential*)malloc(sizeof(Potential));
	
	if(potential)
		fprintf(stdout, "Potential size : %d\n",sizeof(Potential) );
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