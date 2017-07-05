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
				//if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
				//	subcost = 0;
				//}
				//else
				//	subcost = 1;
				subcost = !(((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON);
				int state_ismatch = state[0][threadID] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID+1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_up;
				//去除if的表达方式，是否可以提升性能？
				//myState = (state_ismatch < state_up) * state_ismatch + (state_left < state_up) * state_up + (state_left >= state_up) * state_left;
				
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
					myState = state_up;
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
					myState = state_up;
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


//__global__ void testSharedMemory()
//{
//	__shared__ SPoint queryTraS[MAXLENGTH];
//	__shared__ SPoint traData[MAXLENGTH];
//	__shared__ SPoint traData2[MAXLENGTH];
//	SPoint s;
//	s.x = 4;
//	s.y = 5;
//	traData[1535] = s;
//	queryTraS[1535] = s;
//	traData2[1535] = s;
//}


/*
SPoint版本
同时处理若干个query的EDR，这里以一个EDR计算为单位，每个block计算一个EDR，thread负责一条斜线上state的并行计算。
case1：轨迹长度可长于512，利用循环处理多余512的。
并行计算n个DP
需要提前给定前两次dp的结果，保存在共享内存里
queryTaskNum:总共有几个EDR计算任务
queryTaskInfo[]：每个task对应的qID、candidateID信息，用struct存储
queryTra[],candidateTra[]:轨迹数据，candidateTra保证其内部轨迹不重复
queryTraOffset[],candidateTraOffset[]:每条轨迹的offset，candidateTra保证其内部轨迹不重复
queryLength[],candidateLength[]:每条轨迹的长度（其实offset相减就是长度），其idx和上面的对应
即：candidateLength[id]是第id个candidate Traj的长度
stateTableGPU[][]:对每个candidate的state表
result[]:对于每个candidate的EDR结果
优化方向：
1、轨迹存在share memory里面
2、直接传递轨迹，不使用指针
*/



__global__ void EDRDistance_Batch(int queryTaskNum, TaskInfoTableForSimilarity* taskInfoTable, SPoint *queryTra, int* queryTraOffset, SPoint** candidateTraOffsets, int* queryLength, int *candidateLength, int *result) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	if (blockID >= queryTaskNum) return;
	int thisQueryID = taskInfoTable[blockID].qID;
	int thisQueryLength = queryLength[thisQueryID];
	if ((threadID >= candidateLength[blockID]) && (threadID >= thisQueryLength)) return;
	const int lenT = candidateLength[blockID];
	//int iterNum = queryLength;
	//if (lenT > queryLength)
	//	iterNum = lenT;
	const int iterNum = thisQueryLength + lenT - 1;
	__shared__ int state[2][MAXLENGTH]; //用于存储前两次的结果（占用8KB）
	state[0][0] = 0;
	state[1][0] = 1;
	state[1][1] = 1;
	//对两个轨迹排序，保证第一个比第二个短
	//首先把轨迹存在共享内存里
	//这里面临着share memory是否够用的问题，书上写的是64KB，然而K80似乎有512KB
	//如果是64KB的话，每条轨迹最长1024个点（两个轨迹共占用24KB）
	__shared__ SPoint queryTraS[MAXLENGTH];
	__shared__ SPoint traData[MAXLENGTH];

	for (int i = 0; i <= lenT - 1;i+=MAXTHREAD)
	{
		if(threadID+i<lenT)
		{
			traData[threadID + i] = SPoint(candidateTraOffsets[blockID][threadID + i]);
		}
	}

	SPoint* queryTraBaseAddr = queryTra + queryTraOffset[thisQueryID];
	for (int i = 0; i <= thisQueryLength - 1;i+=MAXTHREAD)
	{
		if(threadID+i<thisQueryLength)
		{
			queryTraS[threadID + i] = *(queryTraBaseAddr + threadID + i);
		}
	}

	const SPoint *tra1, *tra2; //保证tra1比tra2短
	int len1, len2;
	if (lenT >= thisQueryLength) {
		tra1 = queryTraS;
		tra2 = traData;
		len1 = thisQueryLength;
		len2 = lenT;
	}
	else
	{
		tra1 = traData;
		tra2 = queryTraS;
		len1 = lenT;
		len2 = thisQueryLength;
	}

	int myState[5];
	int nodeID;
	for (int i = 0; i <= iterNum - 1; i++) {//第i轮dp
		if (i < len1 - 1) {
			for (int startIdx = 0; startIdx <= i; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID <= i) {
					SPoint p1 = tra1[nodeID];
					SPoint p2 = tra2[i - nodeID]; //这样做内存是聚集访问的吗？
					bool subcost;
					if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
						subcost = 0;
					}
					else
						subcost = 1;
					//subcost = !(((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON);
					int state_ismatch = state[0][nodeID] + subcost;
					int state_up = state[1][nodeID] + 1;
					int state_left = state[1][nodeID + 1] + 1;
					bool c1 = ((state_ismatch < state_up) && (state_ismatch < state_left));
					bool c2 = ((state_left < state_up) && ((state_left < state_ismatch)));
					//去除if的表达方式，是否可以提升性能？
					myState[nodeID / MAXTHREAD] = c1 * state_ismatch + c2 * state_left + !(c1 || c2) * state_up;
					//if ((state_ismatch < state_up) && (state_ismatch < state_left))
					//	myState[nodeID/MAXTHREAD] = state_ismatch;
					//else if ((state_left < state_up) && ((state_left < state_ismatch)))
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					////去除if的表达方式，是否可以提升性能？
					//myState[nodeID / MAXTHREAD] = (state_ismatch < state_up) && (state_ismatch < state_left) * state_ismatch + ((state_left < state_up) && ((state_left < state_ismatch))) * state_left + !(((state_ismatch < state_up) && (state_ismatch < state_left))||(((state_left < state_up) && ((state_left < state_ismatch))))) * state_up;
				}
			}
		}
		else if (i > iterNum - len1) {
			for (int startIdx = 0; startIdx <= iterNum - i - 1; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID <= iterNum - i - 1) {
					SPoint p1 = tra1[nodeID + len1 - (iterNum - i)];
					SPoint p2 = tra2[len2 - 1 - nodeID]; //这样做内存是聚集访问的吗？
					bool subcost;
					if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
						subcost = 0;
					}
					else
						subcost = 1;
					int state_ismatch = state[0][nodeID + 1] + subcost;
					int state_up = state[1][nodeID] + 1;
					int state_left = state[1][nodeID + 1] + 1;
					//if (state_ismatch < state_up)
					//	myState[nodeID / MAXTHREAD] = state_ismatch;
					//else if (state_left < state_up)
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					bool c1 = ((state_ismatch < state_up) && (state_ismatch < state_left));
					bool c2 = ((state_left < state_up) && ((state_left < state_ismatch)));
					//去除if的表达方式，是否可以提升性能？
					myState[nodeID / MAXTHREAD] = c1 * state_ismatch + c2 * state_left + !(c1 || c2) * state_up;
				}
			}
		}
		else
		{
			for (int startIdx = 0; startIdx <= len1; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID <= len1) {
					SPoint p1 = tra1[nodeID];
					SPoint p2 = tra2[i - nodeID]; //这样做内存是聚集访问的吗？
					bool subcost;
					if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
						subcost = 0;
					}
					else
						subcost = 1;
					int state_ismatch = state[0][nodeID] + subcost;
					int state_up = state[1][nodeID] + 1;
					int state_left = state[1][nodeID + 1] + 1;
					//if (state_ismatch < state_up)
					//	myState[nodeID / MAXTHREAD] = state_ismatch;
					//else if (state_left < state_up)
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					bool c1 = ((state_ismatch < state_up) && (state_ismatch < state_left));
					bool c2 = ((state_left < state_up) && ((state_left < state_ismatch)));
					//去除if的表达方式，是否可以提升性能？
					myState[nodeID / MAXTHREAD] = c1 * state_ismatch + c2 * state_left + !(c1 || c2) * state_up;
				}
			}
		}
		//写myState到share内存,ckecked
		
		//int startidx;
		////首先将老数据写到全局内存，全写
		////startidx是旧的数据应该在全局内存中地址，以i-2计算
		////计算应写入全局内存的起始位置
		//// 7.2 发现：似乎stateTableGPU都不用写……是这样吗？
		////////

		//if (i - 2 < len1 - 2) {
		//	startidx = (i - 2 + 2)*(i - 2 + 3) / 2;
		//	for (int Idx = 0; Idx <= i; Idx += MAXTHREAD) {
		//		//if (threadID <= i) {
		//		if(Idx + threadID <= i){
		//			//stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
		//			stateTableGPU[blockID][Idx + threadID + startidx] = state[0][threadID + Idx];
		//		}
		//	}
		//}
		//else if (i - 2 >= iterNum - len1) {
		//	startidx = (len1 + 1)*(len2 + 1) - (iterNum - (i - 2))*(iterNum - (i - 2) + 1) / 2;
		//	for (int Idx = 0; Idx <= iterNum - i + 1; Idx += MAXTHREAD) {
		//		//if (threadID <= iterNum - i + 1) {
		//		if (threadID + Idx <= iterNum - i + 1) {
		//			//stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
		//			stateTableGPU[blockID][Idx + threadID + startidx] = state[0][Idx + threadID];
		//		}
		//	}
		//}
		//else
		//{
		//	startidx = (len1 + 1)*((i - 2) - (len1 - 2)) + len1*(len1 + 1) / 2;
		//	for (int Idx = 0; Idx <= len1; Idx += MAXTHREAD) {
		//		//if (threadID <= len1) {
		//		if (threadID + Idx <= len1) {
		//			// stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
		//			stateTableGPU[blockID][Idx + threadID + startidx] = state[0][Idx + threadID];
		//		}
		//	}
		//}

		//移动新数据到旧数据
		for (int Idx = 0; Idx < MAXLENGTH;Idx+=MAXTHREAD)
		{
			state[0][threadID+Idx] = state[1][threadID+Idx];
		}
		//state[0][threadID] = state[1][threadID];

		//写入新数据
		if (i < len1 - 1) {
			//if (threadID <= i)
			//	state[1][threadID + 1] = myState;
			//if (threadID == 0) {
			//	state[1][0] = i + 2;
			//	state[1][i + 2] = i + 2;
			//}
			for (int Idx = 0; Idx <= i; Idx += MAXTHREAD) {
				if (threadID + Idx <= i)
					state[1][Idx + threadID + 1] = myState[Idx/MAXTHREAD];
			}
			if (threadID == 0) {
				state[1][0] = i + 2;
				state[1][i + 2] = i + 2;
			}
		}
		else if (i >= iterNum - len1) {
			//if (threadID <= iterNum - i - 1)
			//	state[1][threadID] = myState;
			for (int Idx = 0; Idx <= iterNum - i - 1; Idx += MAXTHREAD) {
				if (threadID + Idx <= iterNum - i - 1)
					state[1][threadID + Idx] = myState[Idx / MAXTHREAD];
			}
		}
		else
		{
			//if (threadID < len1)
			//	state[1][threadID + 1] = myState;
			//if (threadID == 0) {
			//	state[1][0] = i + 2;
			//}
			for (int Idx = 0; Idx <= len1; Idx += MAXTHREAD) {
				if (threadID + Idx < len1)
					state[1][Idx + threadID + 1] = myState[Idx / MAXTHREAD];
			}
			if (threadID == 0) {
				state[1][0] = i + 2;
			}
		}
		__syncthreads();
	}
	//输出结果，最后一次计算一定是由进程0完成的
	if (threadID == 0 && blockID < queryTaskNum)
		result[blockID] = myState[0];
}

int EDRDistance_Batch_Handler(int queryTaskNum, TaskInfoTableForSimilarity* taskInfoTable, SPoint *queryTra, int* queryTraOffset, SPoint** candidateTraOffsets, int* queryLength, int *candidateLength, int *result, cudaStream_t *stream)
{
	EDRDistance_Batch <<<queryTaskNum, MAXTHREAD,0 , *stream >>>(queryTaskNum, taskInfoTable, queryTra, queryTraOffset, candidateTraOffsets, queryLength, candidateLength, result);
	return 0;
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
	//
	time1.stop();
	std::cout << time1.elapse() << std::endl;
	time1.start();
	//
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

