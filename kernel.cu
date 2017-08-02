
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cudaKernel.h"
#include "thrust/device_ptr.h"
#include "thrust/remove.h"
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include "ConstDefine.h"


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
__global__ void EDRDistance_1(SPoint *queryTra, SPoint **candidateTra, int candidateNum, int queryLength, int *candidateLength, int** stateTableGPU, int *result) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	if (blockID >= candidateNum) return;
	if ((threadID >= candidateLength[blockID]) && (threadID >= queryLength)) return;
	const int lenT = candidateLength[blockID];
	//int iterNum = queryLength;
	//if (lenT > queryLength)
	//	iterNum = lenT;
	const int iterNum = queryLength + lenT - 1;
	__shared__ short state[2][MAXTHREAD]; //用于存储前两次的结果
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
				int state_left = state[1][threadID + 1] + 1;
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
				SPoint p1 = tra1[threadID + len1 - (iterNum - i)];
				SPoint p2 = tra2[len2 - 1 - threadID]; //这样做内存是聚集访问的吗？
				bool subcost;
				if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
					subcost = 0;
				}
				else
					subcost = 1;
				int state_ismatch = state[0][threadID + 1] + subcost;
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
				SPoint p2 = tra2[i - threadID]; //这样做内存是聚集访问的吗？
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

		if (i - 2 < len1 - 2) {
			startidx = (i - 2 + 2)*(i - 2 + 3) / 2;
			if (threadID <= i) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}
		else if (i - 2 >= iterNum - len1) {
			startidx = (len1 + 1)*(len2 + 1) - (iterNum - (i - 2))*(iterNum - (i - 2) + 1) / 2;
			if (threadID <= iterNum - i + 1) {
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
		if (i < len1 - 1) {
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
	__shared__ int thisQueryID;
	__shared__ int thisQueryLength;
	__shared__ int lenT;
	__shared__ int iterNum;
	//用于存储前两次的结果（占用8KB）

	//对两个轨迹排序，保证第一个比第二个短
	//首先把轨迹存在共享内存里
	//这里面临着share memory是否够用的问题，书上写的是64KB，然而K80似乎有512KB
	//如果是64KB的话，每条轨迹最长1024个点（两个轨迹共占用24KB）
	__shared__ short state[2][MAXLENGTH + 1];
	__shared__ SPoint *queryTraS;
	__shared__ SPoint *traData;
	__shared__ SPoint *tra1, *tra2; //保证tra1比tra2短
	__shared__ int len1, len2;

	if (threadID == 0) {
		thisQueryID = taskInfoTable[blockID].qID;
		thisQueryLength = queryLength[thisQueryID];
		lenT = candidateLength[blockID];
		iterNum = thisQueryLength + lenT - 1;
		state[0][0] = 0;
		state[1][0] = 1;
		state[1][1] = 1;
		queryTraS = queryTra + queryTraOffset[thisQueryID];
		traData = candidateTraOffsets[blockID];
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
	}
	__syncthreads();

	if ((threadID >= lenT) && (threadID >= thisQueryLength)) return;
	



	//__shared__ SPoint queryTraS[MAXLENGTH];
	//__shared__ SPoint traData[MAXLENGTH];


	//for (int i = 0; i <= lenT - 1;i+=MAXTHREAD)
	//{
	//	if(threadID+i<lenT)
	//	{
	//		traData[threadID + i] = SPoint(candidateTraOffsets[blockID][threadID + i]);
	//	}
	//}

	//SPoint* queryTraBaseAddr = queryTra + queryTraOffset[thisQueryID];
	//for (int i = 0; i <= thisQueryLength - 1;i+=MAXTHREAD)
	//{
	//	if(threadID+i<thisQueryLength)
	//	{
	//		queryTraS[threadID + i] = *(queryTraBaseAddr + threadID + i);
	//	}
	//}
	
	int myState[5];
	int nodeID;
	SPoint p1;
	SPoint p2;
	bool subcost;

	for (int i = 0; i <= iterNum - 1; i++) {//第i轮dp
		if (i < len1 - 1) {
			for (int startIdx = 0; startIdx <= i; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID <= i) {
					p1 = tra1[nodeID];
					p2 = tra2[i - nodeID]; //这样做内存是聚集访问的吗？
					if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
						subcost = 0;
					}
					else
						subcost = 1;
					//subcost = !(((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON);
					bool c1 = ((state[0][nodeID] + subcost < (state[1][nodeID] + 1)) && (state[0][nodeID] + subcost < (state[1][nodeID + 1] + 1)));
					bool c2 = (((state[1][nodeID + 1] + 1) < (state[1][nodeID] + 1)) && (((state[1][nodeID + 1] + 1) < state[0][nodeID] + subcost)));
					//去除if的表达方式，是否可以提升性能？
					myState[nodeID / MAXTHREAD] = c1 * (state[0][nodeID] + subcost) + c2 * (state[1][nodeID + 1] + 1) + !(c1 || c2) * (state[1][nodeID] + 1);
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
					p1 = tra1[nodeID + len1 - (iterNum - i)];
					p2 = tra2[len2 - 1 - nodeID]; //这样做内存是聚集访问的吗？
					if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
						subcost = 0;
					}
					else
						subcost = 1;
					//if (state_ismatch < state_up)
					//	myState[nodeID / MAXTHREAD] = state_ismatch;
					//else if (state_left < state_up)
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					bool c1 = (((state[0][nodeID + 1] + subcost) < (state[1][nodeID] + 1)) && ((state[0][nodeID + 1] + subcost) < (state[1][nodeID + 1] + 1)));
					bool c2 = (((state[1][nodeID + 1] + 1) < (state[1][nodeID] + 1)) && (((state[1][nodeID + 1] + 1) < (state[0][nodeID + 1] + subcost))));
					//去除if的表达方式，是否可以提升性能？
					myState[nodeID / MAXTHREAD] = c1 * (state[0][nodeID + 1] + subcost) + c2 * (state[1][nodeID + 1] + 1) + !(c1 || c2) * (state[1][nodeID] + 1);
				}
			}
		}
		else
		{
			for (int startIdx = 0; startIdx < len1; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID < len1) {
					p1 = tra1[nodeID];
					p2 = tra2[i - nodeID]; //这样做内存是聚集访问的吗？
					if (((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)) < EPSILON) {
						subcost = 0;
					}
					else
						subcost = 1;
					//int state_ismatch = (state[0][nodeID] + subcost);
					//int state_up = (state[1][nodeID] + 1);
					//int state_left = (state[1][nodeID + 1] + 1);
					//if (state_ismatch < state_up)
					//	myState[nodeID / MAXTHREAD] = state_ismatch;
					//else if (state_left < state_up)
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					bool c1 = (((state[0][nodeID] + subcost) < (state[1][nodeID] + 1)) && ((state[0][nodeID] + subcost) < (state[1][nodeID + 1] + 1)));
					bool c2 = (((state[1][nodeID + 1] + 1) < (state[1][nodeID] + 1)) && (((state[1][nodeID + 1] + 1) < (state[0][nodeID] + subcost))));
					//去除if的表达方式，是否可以提升性能？
					myState[nodeID / MAXTHREAD] = c1 * (state[0][nodeID] + subcost) + c2 * (state[1][nodeID + 1] + 1) + !(c1 || c2) * (state[1][nodeID] + 1);
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
		for (int Idx = 0; Idx < MAXLENGTH; Idx += MAXTHREAD)
		{
			if(threadID + Idx < MAXLENGTH)
				state[0][threadID + Idx] = state[1][threadID + Idx];
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
					state[1][Idx + threadID + 1] = myState[Idx / MAXTHREAD];
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
	//printf("run kernel now\n");
	EDRDistance_Batch << <queryTaskNum, MAXTHREAD, 0, *stream >> >(queryTaskNum, taskInfoTable, queryTra, queryTraOffset, candidateTraOffsets, queryLength, candidateLength, result);
	return 0;
}

__device__ inline int binary_search_intPair(intPair* temp, int left, int right, int val)
{
	int mid = (left + right) / 2;
	while (left <= right)
	{
		mid = (left + right) / 2;
		if (temp[mid].int_1 == val)
			return temp[mid].int_2;
		else if (temp[mid].int_1 > val)
		{
			right = mid - 1;
		}
		else
			left = mid + 1;
	}
	return 0;
}

__device__ inline int binary_search_intPair_Neighbor(intPair* temp, int left, int right, int val)
{
	int mid = (left + right) / 2;
	while (left <= right)
	{
		mid = (left + right) / 2;
		if (temp[mid].int_1 == val)
			return mid;
		else if (temp[mid].int_1 > val)
		{
			right = mid - 1;
		}
		else
			left = mid + 1;
	}
	return -1;
}

// -1为没找到
__device__ inline int binary_search_int(int* temp, int left, int right, int val)
{
	int mid = (left + right) / 2;
	while (left <= right)
	{
		mid = (left + right) / 2;
		if (temp[mid] == val)
			return mid;
		else if (temp[mid] > val)
		{
			right = mid - 1;
		}
		else
			left = mid + 1;
	}
	return -1;
}

__device__ inline int getIdxFromXYGPU(int x, int y)
{
	int lenx, leny;
	if (x == 0)
		lenx = 1;
	else
	{
		lenx = int(log2f(x)) + 1;
	}
	if (y == 0)
		leny = 1;
	else
		leny = int(log2f(y)) + 1;
	int result = 0;
	int xbit = 1, ybit = 1;
	for (int i = 1; i <= 2 * max(lenx, leny); i++)
	{
		if ((i & 1) == 1) //奇数
		{
			result += (x >> (xbit - 1) & 1) * (1 << (i - 1));
			xbit = xbit + 1;
		}
		else //偶数
		{
			result += (y >> (ybit - 1) & 1) * (1 << (i - 1));
			ybit = ybit + 1;
		}
	}
	return result;
}

__device__ inline int findNeighborGPU(int cellNum, int cellID, int * neighborID)
{
	int x = 0, y = 0;
	for (int bit = 0; bit <= int(log2f(cellNum)) - 1; bit++) {
		if (bit % 2 == 0) {
			//奇数位
			x += ((cellID >> bit)&(1))*(1 << (bit / 2));
		}
		else {
			//偶数位
			y += ((cellID >> bit)&(1))*(1 << (bit / 2));
		}
	}
	int cnt = 0;
	for (int xx = x - 1; xx <= x + 1; xx++) {
		for (int yy = y - 1; yy <= y + 1; yy++) {
			if ((xx != x) || (yy != y))
				neighborID[cnt++] = getIdxFromXYGPU(xx, yy);
			//printf("%d\t", cnt);
		}
	}
	return 0;
}

__device__ inline bool isPositive(short x)
{
	return x >= 0;
}

__global__ void Calculate_FD_Sparse(intPair* queryFVGPU, intPair* FVinfo, intPair* FVTable, intPair* SubbedArray, intPair* SubbedArrayOffset, int SubbedArrayJump, int queryCellLength, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance)
{
	//第一阶段：并行减法
	const int MAX_QUERY_CELLNUMBER = 512;
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int threadIDGlobal = blockDim.x*blockID + threadID;

	__shared__ intPair queryCellTraj[MAX_QUERY_CELLNUMBER];
	__shared__ intPair dbCellTraj[MAX_QUERY_CELLNUMBER];
	//cellchecked记录在query中出现的cell编号，用于在反向减法的时候检查是不是已经减过了。以后可以在归并加中复用此变量。
	__shared__ int cellChecked[MAX_QUERY_CELLNUMBER];
	for (int i = 0; i <= queryCellLength - 1; i += MAXTHREAD) {
		if (threadID + i < queryCellLength)
		{
			queryCellTraj[threadID + i] = queryFVGPU[threadID + i];
		}
	}
	int dbTrajStartIdx = FVinfo[startTrajIdx + blockID].int_2;
	int dbTrajEndIdx;
	if (blockID + startTrajIdx == trajNumInDB - 1)
		dbTrajEndIdx = nonZeroFVNumInDB - 1;
	else
		dbTrajEndIdx = FVinfo[startTrajIdx + blockID + 1].int_2 - 1;

	for (int i = 0; i <= dbTrajEndIdx - dbTrajStartIdx; i += MAXTHREAD)
	{
		if (threadID + i <= dbTrajEndIdx - dbTrajStartIdx)
			dbCellTraj[threadID + i] = FVTable[dbTrajStartIdx + threadID + i];
	}
	//1.1:用query减去db
	for (int i = 0; i < queryCellLength; i += MAXTHREAD)
	{
		if (threadID + i < queryCellLength) {
			int find = binary_search_intPair(dbCellTraj, 0, dbTrajEndIdx - dbTrajStartIdx, queryCellTraj[threadID + i].int_1);
			cellChecked[threadID + i] = queryCellTraj[threadID + i].int_1;
			SubbedArray[SubbedArrayJump * blockID + threadID + i].int_1 = queryCellTraj[threadID + i].int_1;
			SubbedArray[SubbedArrayJump * blockID + threadID + i].int_2 = queryCellTraj[threadID + i].int_2 - find;
		}
		if (threadID == 0) {
			SubbedArrayOffset[blockID].int_1 = queryCellLength - 1;
			SubbedArrayOffset[blockID].int_2 = queryCellLength + dbTrajEndIdx - dbTrajStartIdx;
		}
	}
	//1.2：用db减去query，注意加负号
	for (int i = 0; i <= dbTrajEndIdx - dbTrajStartIdx; i += MAXTHREAD)
	{
		if (threadID + i <= dbTrajEndIdx - dbTrajStartIdx)
		{
			intPair cellNo = dbCellTraj[threadID + i];
			int find = binary_search_int(cellChecked, 0, queryCellLength - 1, cellNo.int_1);
			if (find == -1)
			{
				SubbedArray[SubbedArrayJump * blockID + queryCellLength + threadID + i].int_1 = cellNo.int_1;
				SubbedArray[SubbedArrayJump * blockID + queryCellLength + threadID + i].int_2 = -cellNo.int_2;
			}
			else
				SubbedArray[SubbedArrayJump * blockID + queryCellLength + threadID + i].int_1 = -1;
		}
	}
	__syncthreads();
	//第二阶段：查找相邻，做减法
	//这个阶段改为每个thread处理一个FD
	//2.1：合并每个subbedArray
	if (threadIDGlobal < checkNum) {
		int startMergeIdx = SubbedArrayOffset[threadIDGlobal].int_1 + 1;
		int endMergeIdx = SubbedArrayOffset[threadIDGlobal].int_2;
		int frontPtr = startMergeIdx;
		for (int i = startMergeIdx; i <= endMergeIdx; i++)
		{
			if (SubbedArray[SubbedArrayJump * threadIDGlobal + i].int_1 != -1)
			{
				SubbedArray[SubbedArrayJump * threadIDGlobal + frontPtr] = SubbedArray[SubbedArrayJump * threadIDGlobal + i];
				frontPtr++;
			}
		}
		SubbedArrayOffset[threadIDGlobal].int_2 = frontPtr - 1;
	}
	//2.2 查找相邻
	int neighborsID[8];
	//cell单纯指第几个元素
	for (int cell = 0; cell <= SubbedArrayOffset[threadIDGlobal].int_2; cell++)
	{
		findNeighborGPU(cellNum, cell, neighborsID);
		//for (int i = 0; i <= 7; i++)
		//	neighborsID[i] = 11;
		for (int i = 0; i <= 7; i++)
		{
			int find = binary_search_intPair_Neighbor(&SubbedArray[SubbedArrayJump * threadIDGlobal], 0, SubbedArrayOffset[threadIDGlobal].int_1, neighborsID[i]);
			if (find == -1) {
				find = binary_search_intPair_Neighbor(&SubbedArray[SubbedArrayJump * threadIDGlobal], SubbedArrayOffset[threadIDGlobal].int_1 + 1, SubbedArrayOffset[threadIDGlobal].int_2, neighborsID[i]);
			}
			// 如果是-1，说明这个neighbor是0，不用处理
			if (find != -1)
			{
				if (isPositive(SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2) != isPositive(SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2))
				{
					if (fabsf(SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2) > fabsf(SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2))
					{
						SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2 = SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2 + SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2;
						SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2 = 0;
					}
					else {
						SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2 = SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2 + SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2;
						SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2 = 0;
						break;
					}
				}
			}
		}
	}
	__syncthreads();
	//第三阶段：统计正负个数
	//依然是每个block负责一个FD的计算
	if (blockID >= checkNum)
		return;
	int *tempsumPosi = (int*)queryCellTraj;
	int *tempsumNega = (int*)dbCellTraj;
	tempsumPosi[threadID] = 0;
	tempsumNega[threadID] = 0;
	for (int i = 0; i <= SubbedArrayOffset[blockID].int_2; i += MAXTHREAD)
	{
		if (i + threadID <= SubbedArrayOffset[blockID].int_2)
		{
			tempsumPosi[threadID] += (isPositive(SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2)*SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2);
			tempsumNega[threadID] += (-(!isPositive(SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2))*SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2);
		}
	}
	__shared__ int sizeOfTempSum;
	if (threadID == 0)
		sizeOfTempSum = MAXTHREAD;
	__syncthreads();
	while ((sizeOfTempSum>1))
	{
		if (threadID <= (sizeOfTempSum >> 1) - 1)
		{
			tempsumPosi[threadID] = tempsumPosi[threadID] + tempsumPosi[threadID + (sizeOfTempSum >> 1)];
			tempsumNega[threadID] = tempsumNega[threadID] + tempsumNega[threadID + (sizeOfTempSum >> 1)];
		}
		__syncthreads();
		if (threadID == 0)
			sizeOfTempSum = (sizeOfTempSum >> 1);
		__syncthreads();
	}
	if (threadID == 0)
		FDistance[blockID] = (tempsumPosi[0] > tempsumNega[0]) ? tempsumPosi[0] : tempsumNega[0];
}

//每个block负责一个FD的计算
__global__ void Calculate_FD_NonColumn(short* queryFVGPU, intPair* FVinfo, intPair* FVTable, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance)
{
	//第一阶段：并行减法
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int threadIDGlobal = blockDim.x*blockID + threadID;
	if (blockID >= checkNum)
		return;
	__shared__ intPair taskInfo;
	if (threadID == 0)
		taskInfo = FVinfo[blockID + startTrajIdx];
	int nextCnt;
	if (blockID + startTrajIdx == trajNumInDB - 1)
		nextCnt = nonZeroFVNumInDB;
	else
		nextCnt = FVinfo[blockID + startTrajIdx + 1].int_2;
	__syncthreads();
	for (int i = 0; i <= (cellNum - 1); i += MAXTHREAD)
	{
		int find = binary_search_intPair(FVTable, taskInfo.int_2, (nextCnt - 1), (i + threadID));
		//int find = 1;
		//int k = cellNum*blockID + (i + threadID);
		//queryFVGPU[cellNum*blockID + (i + threadID)] = 2;
		queryFVGPU[cellNum*blockID + (i + threadID)] = queryFVGPU[cellNum*blockID + (i + threadID)] - find;
	}
	//第二阶段：查找相邻，做减法
	//这个阶段改为每个thread处理一个FD
	int neighborsID[8];
	for (int cell = 0; cell <= cellNum - 1; cell++)
	{
		//只需要一部分线程就行了
		if (threadIDGlobal >= checkNum)
			break;
		if (queryFVGPU[cellNum*threadIDGlobal + cell] != 0)
		{
			findNeighborGPU(cellNum, cell, neighborsID);
			//for (int i = 0; i <= 7; i++)
			//	neighborsID[i] = 11;
			for (int i = 0; i <= 7; i++)
			{
				if (isPositive(queryFVGPU[cellNum*threadIDGlobal + cell]) != isPositive(queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]])) {
					if (fabsf(queryFVGPU[cellNum*threadIDGlobal + cell]) > fabsf(queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]]))
					{
						queryFVGPU[cellNum*threadIDGlobal + cell] = queryFVGPU[cellNum*threadIDGlobal + cell] + queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]];
						queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]] = 0;
					}
					else
					{
						queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]] = queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]] + queryFVGPU[cellNum*threadIDGlobal + cell];
						queryFVGPU[cellNum*threadIDGlobal + cell] = 0;
						break;
					}
				}
			}
		}
	}
	__syncthreads();
	//第三阶段：统计正负个数
	//依然是每个block负责一个FD的计算
	__shared__ int tempsumPosi[MAXTHREAD], tempsumNega[MAXTHREAD];
	tempsumPosi[threadID] = 0;
	tempsumNega[threadID] = 0;
	for (int i = 0; i <= cellNum - 1; i += MAXTHREAD)
	{
		tempsumPosi[threadID] += (isPositive(queryFVGPU[blockID*cellNum + (i + threadID)])*queryFVGPU[blockID*cellNum + (i + threadID)]);
		tempsumNega[threadID] += (-(!isPositive(queryFVGPU[blockID*cellNum + (i + threadID)]))*queryFVGPU[blockID*cellNum + (i + threadID)]);
	}
	__shared__ int sizeOfTempSum;
	if (threadID == 0)
		sizeOfTempSum = MAXTHREAD;
	__syncthreads();
	while ((sizeOfTempSum>1))
	{
		if (threadID <= (sizeOfTempSum >> 1) - 1)
		{
			tempsumPosi[threadID] = tempsumPosi[threadID] + tempsumPosi[threadID + (sizeOfTempSum >> 1)];
			tempsumNega[threadID] = tempsumNega[threadID] + tempsumNega[threadID + (sizeOfTempSum >> 1)];
		}
		__syncthreads();
		if (threadID == 0)
			sizeOfTempSum = (sizeOfTempSum >> 1);
		__syncthreads();
	}
	if (threadID == 0)
		FDistance[blockID] = (tempsumPosi[0] > tempsumNega[0]) ? tempsumPosi[0] : tempsumNega[0];

}

//SubbedArrayJump是SubbedArray中每一行有多少个元素，供计算idx用
int Similarity_Pruning_Handler(intPair* queryFVGPU, intPair* FVinfo, intPair* FVTable, intPair* SubbedArray, intPair* SubbedArrayOffset, int SubbedArrayJump, int queryCellLength, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance, cudaStream_t stream)
{
#ifdef NOT_COLUMN_ORIENTED
	Calculate_FD_NonColumn << <checkNum, MAXTHREAD, 0, stream >> >(queryFVGPU, FVinfo, FVTable, startTrajIdx, checkNum, cellNum, trajNumInDB, nonZeroFVNumInDB, FDistance);
#else
	Calculate_FD_Sparse << <checkNum, MAXTHREAD, 0, stream >> >(queryFVGPU, FVinfo, FVTable, SubbedArray, SubbedArrayOffset, SubbedArrayJump, queryCellLength, startTrajIdx, checkNum, cellNum, trajNumInDB, nonZeroFVNumInDB, FDistance);
#endif
	return 0;
}


/*
//先按照能否用一个SM执行一个DP来划分任务，再分别调用两种kernel
//constructing...
可优化：
1、queryTra、queryLength甚至candidateLength可以通过传值的方式直接传递到SM的寄存器，减少全局内存的使用

*/
int handleEDRdistance(SPoint *queryTra, SPoint **candidateTra, int candidateNum, int queryLength, int *candidateLength, int *result) {
	MyTimer time1;
	time1.start();

	int** stateTableGPU = NULL;
	//在GPU中为状态表分配内存
	int** temp = NULL;
	temp = (int**)malloc(sizeof(int*)*candidateNum);
	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMalloc((void**)&temp[i], sizeof(int)*(candidateLength[i] + 1)*(queryLength + 1)));
	}
	CUDA_CALL(cudaMalloc((void***)&stateTableGPU, sizeof(int*)*candidateNum));
	CUDA_CALL(cudaMemcpy(stateTableGPU, temp, candidateNum*sizeof(int*), cudaMemcpyHostToDevice));

	//为存储的轨迹信息分配内存
	SPoint *queryTraGPU = NULL, **candidateTraGPU = NULL;
	int *candidateLengthGPU = NULL, *resultGPU = NULL;
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
	EDRDistance_1 << <candidateNum, MAXTHREAD >> >(queryTraGPU, candidateTraGPU, candidateNum, queryLength, candidateLengthGPU, stateTableGPU, resultGPU);

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


inline void __getLastCudaError(const char *errorMessage, const char *file, const int line)
{
	cudaError_t err = cudaGetLastError();

	if (cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
			file, line, errorMessage, (int)err, cudaGetErrorString(err));
			exit(EXIT_FAILURE);
	}
}


//using namespace thrust;
//static const int MAXTHREAD = 512; //每个block线程数

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
void CUDAwarmUp() {
	CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));
	CUDA_CALL(cudaSetDevice(0));
	
}


#ifdef _CELL_BASED_STORAGE
int putCellDataSetIntoGPU(Point* pointsPtr, Point*& pointsPtrGPU, int pointNum) {
	
	CUDA_CALL(cudaMalloc((void**)&pointsPtrGPU, pointNum * sizeof(Point))); //分配数据的内存
	//debug
	//std::cout << pointNum << std::endl;
	//debug
	CUDA_CALL(cudaMemcpy(pointsPtrGPU, pointsPtr, pointNum * sizeof(Point), cudaMemcpyHostToDevice));//数据拷贝到gpu里
	return 0;
}
__global__ void cudaRangeQuery(int* rangeStarts, int* rangeEnds, int candidateCellNum, const Point* pointsPtr, const float xmin, const float ymin, const float xmax, const float ymax, const int *resultOffset, Point* resultPtrCuda) {
	int cellNo = blockIdx.x; //candidate里面第几个cell 0,1,2,....
	if (cellNo >= candidateCellNum) return;
	int tid = threadIdx.x;
	if (tid >= 256) return;
	int pointNum = rangeEnds[cellNo] - rangeStarts[cellNo] + 1;//block要处理的这个cell有这么多个点
	const int offset = rangeStarts[cellNo];
	for (int i = tid; i <= pointNum - 1; i += MAXTHREAD) {
		float x = pointsPtr[offset + i].x;
		float y = pointsPtr[offset + i].y;
		uint32_t tid = pointsPtr[offset + i].tID;
		uint32_t time = pointsPtr[offset + i].time;
		if (x <= xmax &&x >= xmin&&y <= ymax&&y >= ymin) {
			resultPtrCuda[resultOffset[cellNo] + i].x = x;
			resultPtrCuda[resultOffset[cellNo] + i].y = y;
			resultPtrCuda[resultOffset[cellNo] + i].tID = tid;
			resultPtrCuda[resultOffset[cellNo] + i].time = time;
		}
		else
			resultPtrCuda[resultOffset[cellNo] + i].tID = -1;
	}
}

__global__ void cudaRangeQueryTest(RangeQueryStateTable* stateTable, int stateTableLength, uint8_t* result, 
	const int maxTrajNum) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	__shared__ RangeQueryStateTable sharedStateTable;
	// __shared__ uint8_t resultTemp[10000]; //10K
	if(tID == 0)
		sharedStateTable = (stateTable[bID]);
	__syncthreads();
	int jobID = sharedStateTable.queryID;
	SPoint *baseAddr = (SPoint*)(sharedStateTable.ptr);
	int candidateNum = sharedStateTable.candidatePointNum;//该block的需要查询的点的个数
	//int resultOffset = bID*maxPointNumInStateTable; //该block的结果的起始地址
	SPoint p;
	/*
	for (int i = 0; i <= candidateNum / MAXTHREAD-1; i++) {
		//p = *(baseAddr + (i*MAXTHREAD + tID));
		p = baseAddr[i*MAXTHREAD + tID];
		//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
			//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//如果验证通过，则该值为本身编号，否则为0
		//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
		//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
		//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
		if((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
			result[jobID*maxTrajNum + p.tID] = 1;
		//如果验证通过，则相应位被置为1

		//__syncthreads();
	}
	if (tID < candidateNum - candidateNum / MAXTHREAD * MAXTHREAD) {
		p = *(baseAddr + (candidateNum / MAXTHREAD * MAXTHREAD + tID));
		//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
		//	(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(candidateNum / MAXTHREAD * MAXTHREAD + tID);//如果验证通过，则该值为本身编号，否则为0
		//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].jobID = bID;
		//result[resultOffset + (candidateNum / MAXTHREAD * MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) &&
		//	(p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
		if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
			result[jobID*maxTrajNum + p.tID] = 1;
	}
	*/
	
	//new version
	for (int i = 0; i < candidateNum; i+=MAXTHREAD) {
		//p = *(baseAddr + (i*MAXTHREAD + tID));
		if (i + tID < candidateNum) {
			p = baseAddr[i + tID];
			//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
			//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//如果验证通过，则该值为本身编号，否则为0
			//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
			//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
			//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
			if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
				result[jobID*maxTrajNum + p.tID] = 1;
			//如果验证通过，则相应位被置为1
		}
		//__syncthreads();
	}
	


	//else {
	//	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = 0; //多出来的部分，直接设为无效即可
	//}
	//__syncthreads();
	//__syncthreads();
	//int globalTID = blockDim.x * blockIdx.x + threadIdx.x;
	//if (globalTID < stateTableLength) {

	//}
}

__global__ void cudaRangeQuerySTIG(RangeQueryStateTable* stateTable, int stateTableLength, uint8_t* result,
	const int maxTrajNum) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	//if (bID > stateTableLength)
	//	return;
	__shared__ RangeQueryStateTable sharedStateTable;
	// __shared__ uint8_t resultTemp[10000]; //10K
	if (tID == 0)
		sharedStateTable = (stateTable[bID]);
	__syncthreads();
	int jobID = sharedStateTable.queryID;
	SPoint *baseAddr = (SPoint*)(sharedStateTable.ptr);
	int candidateNum = sharedStateTable.candidatePointNum;//该block的需要查询的点的个数
														  //int resultOffset = bID*maxPointNumInStateTable; //该block的结果的起始地址
	SPoint p;
	//new version
	for (int i = 0; i < candidateNum; i += MAXTHREAD) {
		//p = *(baseAddr + (i*MAXTHREAD + tID));
		if (i + tID < candidateNum) {
			p = baseAddr[i + tID];
			//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
			//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//如果验证通过，则该值为本身编号，否则为0
			//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
			//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
			//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
			if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
				result[jobID*maxTrajNum + p.tID] = 1;
			//如果验证通过，则相应位被置为1
		}
		//__syncthreads();
	}



	//else {
	//	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = 0; //多出来的部分，直接设为无效即可
	//}
	//__syncthreads();
	//__syncthreads();
	//int globalTID = blockDim.x * blockIdx.x + threadIdx.x;
	//if (globalTID < stateTableLength) {

	//}
}

int cudaRangeQuerySTIGHandler(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxQueryNum, cudaStream_t stream)
{
	//RangeQueryResultGPU* resultGPU;
	//MyTimer timer;
	uint8_t* resultGPU;
	//int resultByteNum = (maxPointNum)/8+1;//每个结果需要用几个byte保存，不能按比特，只能按字节
	CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxTrajNum)*maxQueryNum));//selective低一点
	CUDA_CALL(cudaMemset(resultGPU, 0, (maxTrajNum)*maxQueryNum));
	//timer.start();
	//多分配一点内存，每个stateTable项占据的内存数相等
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxPointNum)*stateTableLength));

	//CUDA_CALL(cudaMalloc((void**)&resultGPU, maxPointNum*stateTableLength*sizeof(RangeQueryResultGPU)));
	//timer.stop();
	//std::cout << "Time 1:" << timer.elapse() << "ms" << std::endl;

	//timer.start();	
	cudaRangeQuerySTIG << <stateTableLength, MAXTHREAD, 0, stream >> >(stateTableGPU, stateTableLength, resultGPU, maxTrajNum);
	CUDA_CALL(cudaDeviceSynchronize());
	//timer.stop();
	//std::cout << "Time 2:" << timer.elapse() << "ms" << std::endl;

	//timer.start();

	CUDA_CALL(cudaMemcpy(result, resultGPU, (maxTrajNum)*maxQueryNum, cudaMemcpyDeviceToHost));

	//timer.stop();
	//std::cout << "Time 3:" << timer.elapse() << "ms" << std::endl;
	return 0;
}

int cudaRangeQueryTestHandler(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxJobNum, cudaStream_t stream) {
	//RangeQueryResultGPU* resultGPU;
	//MyTimer timer;
	uint8_t* resultGPU;
	//int resultByteNum = (maxPointNum)/8+1;//每个结果需要用几个byte保存，不能按比特，只能按字节
	CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxTrajNum)*maxJobNum));//selective低一点
	CUDA_CALL(cudaMemset(resultGPU, 0, (maxTrajNum)*maxJobNum));
	//timer.start();
	//多分配一点内存，每个stateTable项占据的内存数相等
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxPointNum)*stateTableLength));
	
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, maxPointNum*stateTableLength*sizeof(RangeQueryResultGPU)));
	//timer.stop();
	//std::cout << "Time 1:" << timer.elapse() << "ms" << std::endl;
	
	//timer.start();	
	cudaRangeQueryTest <<<stateTableLength, MAXTHREAD,0, stream >>>(stateTableGPU, stateTableLength, resultGPU, maxTrajNum);
	CUDA_CALL(cudaDeviceSynchronize());
	//timer.stop();
	//std::cout << "Time 2:" << timer.elapse() << "ms" << std::endl;

	//timer.start();
	
	CUDA_CALL(cudaMemcpy(result, resultGPU, (maxTrajNum)*maxJobNum, cudaMemcpyDeviceToHost));
	CUDA_CALL(cudaFree(resultGPU));

	//timer.stop();
	//std::cout << "Time 3:" << timer.elapse() << "ms" << std::endl;
	return 0;
}

/*
int cudaRangeQueryHandler(int* candidateCells, int* rangeStarts, int* rangeEnds, int candidateCellNum,float xmin, float ymin, float xmax, float ymax, Point*& resultsGPU, int& resultNum,Point *pointsPtrGPU,Point *&result) {
	//第一个参数暂时没有用，注意这里candidatecells[i]已经不再是cell的id，由于只存非空
	//第四个参数表示非空的cell个数
	//计算candidate集中的点的总数，gpu内开辟相同大小的空间做flag。rangestart和rangeend是相应candidatecell内的采样点在AllPoints的起始下标和终止下标
	//倒数第三个和倒数第二个参数是输出，一个是保存结果的GPU地址，第二个是结果的个数
	//PointsPtrGPU是数据集在gpu的地址
	MyTimer timer1;
	timer1.start();
	int counter = 0;
	int *resultOffset = (int*)malloc(sizeof(int)*candidateCellNum);
	//std::cout << candidateCellNum << ":"<<std::endl;
	for (int i = 0; i <= candidateCellNum - 1; i++) {
		resultOffset[i] = counter;
		////debug
		//std::cout << "(" << rangeStarts[i] << "," << rangeEnds[i] << ");"<<"["<<resultOffset[i]<<"]";
		////debug
		counter += rangeEnds[i] - rangeStarts[i] + 1;
	}
	int totalPointNumInCandidate = counter;


	int *rangeStartsCuda = NULL, *rangeEndsCuda = NULL, *resultOffsetCuda = NULL;

	CUDA_CALL(cudaMalloc((void**)&resultsGPU, sizeof(Point)*totalPointNumInCandidate));
	//将range和cell信息写入gpu
	//CUDA_CALL(cudaMalloc((void**)&candidateCellsCuda, sizeof(int)*candidateCellNum));
	CUDA_CALL(cudaMalloc((void**)&rangeStartsCuda, candidateCellNum*sizeof(int)));
	//std::cout << "\n" << candidateCellNum*sizeof(int) << "\n";
	CUDA_CALL(cudaMalloc((void**)&rangeEndsCuda, candidateCellNum*sizeof(int)));
	CUDA_CALL(cudaMalloc((void**)&resultOffsetCuda, candidateCellNum*sizeof(int)));
	//CUDA_CALL(cudaMemcpy(candidateCellsCuda, candidateCells, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(rangeStartsCuda, rangeStarts, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMemcpy(rangeEndsCuda, rangeEnds, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(resultOffsetCuda, resultOffset, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	////debug
	//CUDA_CALL(cudaMemcpy( rangeStarts, rangeStartsCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaMemcpy( rangeEnds, rangeEndsCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaMemcpy( resultOffset, resultOffsetCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//for (int i = 0; i <= candidateCellNum - 1; i++) {
	//	//debug
	//	std::cout << "(" << rangeStarts[i] << "," << rangeEnds[i] << ");" << "[" << resultOffset[i] << "]";
	//	//debug
	//}
	//debug
	timer1.stop();
	std::cout << timer1.ticks() << std::endl;
	timer1.start();
	//调用kernel，如果某点满足条件，在相应位置写入在AllPoints中的下标，否则写入-1
	//每个cell分配给一个block
	cudaRangeQuery <<<candidateCellNum, MAXTHREAD >>>(rangeStartsCuda, rangeEndsCuda, candidateCellNum, pointsPtrGPU, xmin, ymin, xmax, ymax, resultOffsetCuda, resultsGPU);
	//kernel调用结束，结果保存在idxsGPU中，如果符合条件，对应元素内容是在AllPoints的下标，如果不符合，内容为-1
	//CUDA_CALL(cudaFree(candidateCellsCuda));
	CUDA_CALL(cudaFree(rangeStartsCuda));
	CUDA_CALL(cudaFree(rangeEndsCuda));
	CUDA_CALL(cudaFree(resultOffsetCuda));
	//getLastCudaError("Error in Calling 'kernel'");
	//使用Thrust删除所有-1，得到最终结果
	timer1.stop();
	std::cout << timer1.ticks() << std::endl;

	//串行结果合并
	//test

	timer1.start();
	Point *resultset = NULL;
	resultset = (Point*)malloc(totalPointNumInCandidate*sizeof(Point));
	CUDA_CALL(cudaMemcpy(resultset, resultsGPU, sizeof(Point)*totalPointNumInCandidate, cudaMemcpyDeviceToHost));
	std::vector<Point> *resultPoint = new std::vector<Point>;
	for (int i = 0; i <= totalPointNumInCandidate - 1; i++) {
		if (resultset[i].tID != -1)
		{
			resultPoint->push_back(resultset[i]);
		}
	}
	result = &resultPoint->at(0);
	free(resultset);
	//test
	timer1.stop();
	std::cout << timer1.ticks() << std::endl;
	
	//并行结果合并
	
	//thrust::device_ptr<int> idxsPtr = thrust::device_pointer_cast(idxsGPU);
	//int a;
	//cudaMemcpy(&a, idxsGPU, 1, cudaMemcpyDeviceToHost);
	//size_t num = thrust::remove(idxsPtr, idxsPtr + totalPointNumInCandidate-1, -1) - idxsPtr;
	//int *result = (int*)malloc(sizeof(int)*num);
	//thrust::copy(idxsPtr, idxsPtr + num, result);
	//resultNum = num;
	//resultIdx = result;

	//CUDA_CALL(cudaFree(idxsGPU));


	return 0;
}
*/




#else
int cudaRangeQueryHandler(Point* pointsPtr, int pointNum, float xmin, float ymin, float xmax, float ymax, Point*& resultsPtr, int& resultNum) {
	Point* pointsPtrCuda = NULL;
	Point* resultPtrCuda = NULL;
	CUDA_CALL(cudaMalloc((void**)&pointsPtrCuda, pointNum * sizeof(Point))); //分配数据的内存
	CUDA_CALL(cudaMalloc((void**)&resultPtrCuda, pointNum * sizeof(Point))); //gpu内存储结果的地方
	CUDA_CALL(cudaMemcpy(pointsPtrCuda, pointsPtr, pointNum * sizeof(Point), cudaMemcpyHostToDevice));//数据拷贝到gpu里

																									  //调用核函数处理数据，结果放在gpu中

																									  //取回数据，处理
	return 0;
}
#endif




__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

//int main()
//{
//    const int arraySize = 5;
//    const int a[arraySize] = { 1, 2, 3, 4, 5 };
//    const int b[arraySize] = { 10, 20, 30, 40, 50 };
//    int c[arraySize] = { 0 };
//
//    // Add vectors in parallel.
//    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "addWithCuda failed!");
//        return 1;
//    }
//
//    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
//        c[0], c[1], c[2], c[3], c[4]);
//
//    // cudaDeviceReset must be called before exiting in order for profiling and
//    // tracing tools such as Nsight and Visual Profiler to show complete traces.
//    cudaStatus = cudaDeviceReset();
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaDeviceReset failed!");
//        return 1;
//    }
//
//    return 0;
//}

// Helper function for using CUDA to add vectors in parallel.
