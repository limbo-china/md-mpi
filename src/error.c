#include "error.h"

#include <stdio.h>

const char* errInfo[errNums] ={
	"no error",
	"xProcNum * yProcNum * zProcNum != rankNum",

};

// 提示错误信息
void errorInfo(const enum ErrorPtr ptr){

	fprintf(stdout, "error %d: %s\n", ptr, errInfo[ptr]);
}