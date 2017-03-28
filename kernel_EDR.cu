//测试EDR的GPU并行算法
//zbw0046 3.22



#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <stdio.h>
#include "ConstDefine.h"
#include "cudaKernel.h"
#include <assert.h>
#include"device_functions.h"
#include "WinTimer.h"

#define CUDA_CALL(x) { const cudaError_t a = (x); if (a!= cudaSuccess) { printf("\nCUDA Error: %s(err_num=%d)\n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}


/*
并行计算1个大规模dp
需要提前给定前两次dp的结果，保存在共享内存里
iter: 第几个dp单位；outputIdx：输出结果在全局内存位置；tra1、tra2：两条轨迹，提前被载入共享内存；
*/
//__global__ void DPforward(const int iter, const int* outputIdx,const SPoint *tra1,const SPoint *tra2) {
//	SPoint p1 = tra1[threadIdx.x];
//	SPoint p2 = tra2[iter - threadIdx.x - 1]; //这样做内存是聚集访问的吗？
//	bool subcost;
//	if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
//		subcost = 0;
//	}
//	else
//		subcost = 1;
//
//}

/*
SPoint版本
case1：轨迹长度小于512
并行计算n个DP
需要提前给定前两次dp的结果，保存在共享内存里
queryTra[],candidateTra[][]:轨迹
stateTableGPU[][]:对每个candidate的state表
result[]:对于每个candidate的EDR结果
优化方向：
1、轨迹存在share memory里面
2、直接传递轨迹，不使用指针
*/
__global__ void EDRDistance_1(SPoint *queryTra, SPoint **candidateTra,int candidateNum,int queryLength,int *candidateLength,int** stateTableGPU,int *result) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	if (blockID >= candidateNum) return;
	if ((threadID >= candidateLength[blockID]) && (threadID >= queryLength)) return;
	const int lenT = candidateLength[blockID];
	//int iterNum = queryLength;
	//if (lenT > queryLength)
	//	iterNum = lenT;
	const int iterNum = queryLength + lenT - 1;
	__shared__ int state[2][MAXTHREAD]; //用于存储前两次的结果
	state[0][0] = 0;
	state[1][0] = 1;
	state[1][1] = 1;
	//对两个轨迹排序，保证第一个比第二个短
	//首先把轨迹存在共享内存里
	__shared__ SPoint queryTraS[MAXTHREAD];
	__shared__ SPoint traData[MAXTHREAD];
	if (threadID < lenT) {
		traData[threadID] = candidateTra[blockID][threadID];
	}
	if (threadID < queryLength) {
		queryTraS[threadID] = queryTra[threadID];
	}
	const SPoint *tra1, *tra2; //保证tra1比tra2短
	int len1, len2;
	if (lenT >= queryLength) {
		tra1 = queryTraS;
		tra2 = traData;
		len1 = queryLength;
		len2 = lenT;
	}
	else
	{
		tra1 = traData;
		tra2 = queryTraS;
		len1 = lenT;
		len2 = queryLength;
	}

	int myState;
	for (int i = 0; i <= iterNum - 1; i++) {//第i轮dp
		if (i < len1 - 1) {
			if (threadID <= i) {
				SPoint p1 = tra1[threadID];
				SPoint p2 = tra2[i - threadID]; //这样做内存是聚集访问的吗？
				bool subcost;
				if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
					subcost = 0;
				}
				else
					subcost = 1;
				int state_ismatch = state[0][threadID] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID+1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_ismatch;
			}
		}
		else if (i > iterNum - len1) {
			if (threadID <= iterNum - i - 1) {
				SPoint p1 = tra1[threadID+len1-(iterNum-i)];
				SPoint p2 = tra2[len2-1-threadID]; //这样做内存是聚集访问的吗？
				bool subcost;
				if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
					subcost = 0;
				}
				else
					subcost = 1;
				int state_ismatch = state[0][threadID+1] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID + 1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_ismatch;
			}
		}
		else
		{
			if (threadID < len1) {
				SPoint p1 = tra1[threadID];
				SPoint p2 = tra2[i-threadID]; //这样做内存是聚集访问的吗？
				bool subcost;
				if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
					subcost = 0;
				}
				else
					subcost = 1;
				int state_ismatch = state[0][threadID] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID + 1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_ismatch;
			}
		}
		//写myState到share内存,ckecked
		int startidx;
		//首先将老数据写到全局内存，全写
		//startidx是旧的数据应该在全局内存中地址，以i-2计算
		//计算应写入全局内存的起始位置

		if (i-2 < len1 - 2) {
			startidx = (i-2 + 2)*(i-2 + 3) / 2;
			if (threadID <= i) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}
		else if (i-2 >= iterNum - len1) {
			startidx = (len1 + 1)*(len2 + 1) - (iterNum - (i-2))*(iterNum - (i-2) + 1) / 2;
			if (threadID <= iterNum - i + 1 ) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}
		else
		{
			startidx = (len1 + 1)*((i - 2) - (len1 - 2)) + len1*(len1 + 1) / 2;
			if (threadID <= len1) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}

		//移动新数据到旧数据
		state[0][threadID] = state[1][threadID];
		//写入新数据
		if (i < len1-1) {
			if (threadID <= i)
				state[1][threadID + 1] = myState;
			if (threadID == 0) {
				state[1][0] = i + 2;
				state[1][i + 2] = i + 2;
			}
		}
		else if (i >= iterNum - len1) {
			if (threadID <= iterNum - i - 1)
				state[1][threadID] = myState;
		}
		else
		{
			if (threadID < len1)
				state[1][threadID + 1] = myState;
			if (threadID == 0) {
				state[1][0] = i + 2;
			}
		}
		__syncthreads();
	}
	//输出结果，最后一次计算一定是由进程0完成的
	if (threadID == 0)
		result[blockID] = myState;
}
/*
//先按照能否用一个SM执行一个DP来划分任务，再分别调用两种kernel
//constructing...
可优化：
1、queryTra、queryLength甚至candidateLength可以通过传值的方式直接传递到SM的寄存器，减少全局内存的使用

*/
int handleEDRdistance(SPoint *queryTra, SPoint **candidateTra, int candidateNum, int queryLength, int *candidateLength,int *result) {
	MyTimer time1;
	time1.start();

	int** stateTableGPU=NULL;
	//在GPU中为状态表分配内存
	int** temp=NULL;
	temp = (int**)malloc(sizeof(int*)*candidateNum);
	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMalloc((void**)&temp[i], sizeof(int)*(candidateLength[i] + 1)*(queryLength + 1)));
	}
	CUDA_CALL(cudaMalloc((void***)&stateTableGPU, sizeof(int*)*candidateNum));
	CUDA_CALL(cudaMemcpy(stateTableGPU, temp, candidateNum*sizeof(int*), cudaMemcpyHostToDevice));

	//为存储的轨迹信息分配内存
	SPoint *queryTraGPU=NULL, **candidateTraGPU=NULL;
	int *candidateLengthGPU=NULL, *resultGPU=NULL;
	CUDA_CALL(cudaMalloc((void**)&queryTraGPU, sizeof(SPoint)*queryLength));
	CUDA_CALL(cudaMalloc((void**)&candidateLengthGPU, sizeof(int)*candidateNum));
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, sizeof(int)*candidateNum));

	SPoint **tempS = (SPoint**)malloc(sizeof(SPoint*)*candidateNum);
	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMalloc((void**)&tempS[i], sizeof(SPoint)*candidateLength[i]));
		
	}
	CUDA_CALL(cudaMalloc((void***)&candidateTraGPU, sizeof(SPoint*)*candidateNum));
	CUDA_CALL(cudaMemcpy(candidateTraGPU, tempS, candidateNum*sizeof(SPoint*), cudaMemcpyHostToDevice));

	//最好通过传参数的方法传递轨迹，这就要求轨迹连续存储
	//向GPU传递轨迹信息
	CUDA_CALL(cudaMemcpy(queryTraGPU, queryTra, queryLength*sizeof(SPoint), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(candidateLengthGPU, candidateLength, candidateNum*sizeof(int), cudaMemcpyHostToDevice));
	
	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMemcpy(tempS[i], candidateTra[i], candidateLength[i] * sizeof(SPoint), cudaMemcpyHostToDevice));
	}
	//for (int i = 0; i <= candidateNum - 1;i++)
	//	CUDA_CALL(cudaMemcpy(candidateTraGPU[i], candidateTra[i], candidateLength[i]*sizeof(SPoint), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaHostAlloc((void**)&result, candidateNum*sizeof(int), cudaHostAllocWriteCombined | cudaHostAllocMapped));
	CUDA_CALL(cudaHostGetDevicePointer(&resultGPU, result, 0));
	time1.stop();
	std::cout << time1.elapse() << std::endl;
	time1.start();
	//执行kernel
	EDRDistance_1 <<<candidateNum, MAXTHREAD >>>(queryTraGPU, candidateTraGPU, candidateNum, queryLength, candidateLengthGPU, stateTableGPU, resultGPU);

	//取结果
	//result = (int*)malloc(candidateNum*sizeof(int));
	//CUDA_CALL(cudaMemcpy(result, resultGPU, candidateNum*sizeof(int), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();
//	for (int j = 0; j <= candidateNum - 1;j++)
//		std::cout << result[j] << std::endl;

	//free GPU!!!!!
	time1.stop();
	std::cout << time1.elapse() << std::endl;
	return 0;

}

