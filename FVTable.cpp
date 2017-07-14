#include "FVTable.h"
#include <queue>
#include "WinTimer.h"
#include "cudaKernel.h"


int FVTable::initFVTable(int trajNum,int cellNum)
{
	this->cellNum = cellNum;
	this->trajNum = trajNum;
	this->FreqVector.resize(trajNum);
	return 0;
}



//向FVTable中添加pointNum个点
int FVTable::addPointToFVTable(int trajID, int pointNum, int cellID)
{
	map<int, int>::iterator iter = this->FreqVector[trajID].find(cellID);
	if (iter == this->FreqVector[trajID].end()) {
		//没有找到，添加之
		this->FreqVector[trajID].insert(pair<int, int>(cellID, pointNum));
	}
	else {
		//找到，+1
		iter->second = iter->second + pointNum;
	}
	return 0;
}



int FVTable::getCandidate(int bestNow, int k, map<int,int> * freqVectorQ, int * candidateTrajID, int * candidateNum)
{
	priority_queue<FDwithID, vector<FDwithID>, cmp> FDqueue;
	this->formPriorityQueue(&FDqueue, freqVectorQ);
	for (int i = 0; i <= k - 1; i++) {
		candidateTrajID[i] = FDqueue.top().traID;
		if (candidateTrajID[i] >= bestNow)
			return 1;//没有找够k个
		FDqueue.pop();
		(*candidateNum)++;
	}
	return 0;
}

double FVTable::calculateFreqDist(int * freqVectorQ, int trajID)
{
	return 0.0;
}

int FVTable::findNeighbor(int cellID, int * neighborID)
{
	int x = 0, y = 0;
	for (int bit = 0; bit <= int(log2(this->cellNum)) - 1; bit++) {
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
				neighborID[cnt++] = getIdxFromXY(xx, yy);
			//printf("%d\t", cnt);
		}
	}
	return 0;
}

int FVTable::formPriorityQueue(priority_queue<FDwithID, vector<FDwithID>, cmp> *queue, map<int, int>* freqVectorQ)
{
	MyTimer time1;
	for (int i = 1; i <= this->trajNum - 1; i++) {
		//对于该轨迹，计算与查询的FV的FD
		//首先，计算两个vector的对应元素减法
		//原始轨迹减去查询轨迹
		time1.start();
		int tempVector;// = (int*)malloc(sizeof(int)*this->cellNum);
		map<int, int> tempPositive;
		map<int, int> tempNegative;
		for (map<int, int>::iterator iter = this->FreqVector[i].begin(); iter != this->FreqVector[i].end(); iter++) {
			int cid = iter->first;
			int cfreq = iter->second;
			map<int, int>::iterator iter_query = freqVectorQ->find(cid);
			if (iter_query == freqVectorQ->end()) {
				//说明在query里面该频率为0
				tempPositive.insert(pair<int, int>(cid, cfreq));
			}
			else {
				tempVector = cfreq - iter_query->second;
				if (tempVector>0)
					tempPositive.insert(pair<int, int>(cid, tempVector));
				else if (tempVector<0)
					tempNegative.insert(pair<int, int>(cid, -tempVector));
			}
		}
		time1.stop();
		printf("prun time 1:%f\n", time1.elapse());
		time1.start();
		//对于在原始轨迹中没有的
		for (map<int, int>::iterator iter = freqVectorQ->begin(); iter != freqVectorQ->end(); iter++) {
			int cid = iter->first;
			//如果这个cid在上一个步骤处理过了，则跳过
			map<int, int>::iterator iter_to_database = this->FreqVector[i].find(cid);
			if (iter_to_database != this->FreqVector[i].end())
				continue;

			//到这里才说明刚才的计算中没有包含这个cell，但是这个cell在查询轨迹中出现了，原始轨迹中没有
			tempVector = -(iter->second);
			if (tempVector>0)
				tempPositive.insert(pair<int, int>(cid, tempVector));
			else if (tempVector<0)
				tempNegative.insert(pair<int, int>(cid, -tempVector));

		}
		time1.stop();
		printf("prun time 2:%f\n", time1.elapse());
		time1.start();
		//减法完成
		//对正负两个map对邻接的cell处理
		for (map<int, int>::iterator iter = tempPositive.begin(); iter != tempPositive.end(); iter++) {
			//找邻接的cell
			int cid = iter->first;
			int posiValue = iter->second;
			int neighborIDs[8];
			this->findNeighbor(cid, neighborIDs);
			map<int, int>::iterator iter_neigh;
			for (int i = 0; i <= 7; i++) {
				iter_neigh = tempNegative.find(neighborIDs[i]);
				if (iter_neigh != tempNegative.end()) {
					int negaValue = iter_neigh->second;
					if (posiValue > negaValue) {
						tempPositive[cid] = posiValue - negaValue;
						tempNegative.erase(neighborIDs[i]);
					}
					else if (posiValue == negaValue) {
						tempPositive.erase(cid);
						tempNegative.erase(neighborIDs[i]);
						break;
					}
					else {
						tempNegative[neighborIDs[i]] = negaValue - posiValue;
						tempPositive.erase(cid);
						break;
					}
				}
			}
		}
		time1.stop();
		printf("prun time 3:%f\n", time1.elapse());
		time1.start();
		for (map<int, int>::iterator iter = tempNegative.begin(); iter != tempNegative.end(); iter++) {
			//找邻接的cell
			int cid = iter->first;
			int negaValue = iter->second;
			int neighborIDs[8];
			this->findNeighbor(cid, neighborIDs);
			map<int, int>::iterator iter_neigh;
			for (int i = 0; i <= 7; i++) {
				iter_neigh = tempPositive.find(neighborIDs[i]);
				if (iter_neigh != tempPositive.end()) {
					int posiValue = iter_neigh->second;
					if (negaValue > posiValue) {
						tempNegative[cid] = negaValue - posiValue;
						tempPositive.erase(neighborIDs[i]);
					}
					else if (posiValue == negaValue) {
						tempNegative.erase(cid);
						tempPositive.erase(neighborIDs[i]);
						break;
					}
					else {
						tempPositive[neighborIDs[i]] = posiValue - negaValue;
						tempNegative.erase(cid);
						break;
					}
				}
			}
		}
		time1.stop();
		printf("prun time 4:%f\n", time1.elapse());
		time1.start();
		//邻接情况处理完成
		int sumPosi = 0, sumNega = 0;
		for (map<int, int>::iterator iter = tempPositive.begin(); iter != tempPositive.end(); iter++)
			sumPosi += iter->second;
		for (map<int, int>::iterator iter = tempNegative.begin(); iter != tempNegative.end(); iter++)
			sumNega += iter->second;
		//统计正负，取最大值的绝对值
		int resultLB = max(sumPosi, sumNega);
		FDwithID fd;
		fd.traID = i;
		fd.FD = resultLB;
		queue->push(fd);
		time1.stop();
		printf("prun time 5:%f\n", time1.elapse());
	}
	return 0;
}

int FVTable::transferFVtoGPU()
{
	CUDA_CALL(cudaMalloc((void**)&this->FVinfoGPU, 16 * 1024 * 1024));
	CUDA_CALL(cudaMalloc((void**)&this->FVTableOffset, 256 * 1024 * 1024));
	//列式存储，容量和之前offset相同
	CUDA_CALL(cudaMalloc((void**)&this->FVTableGPU, 256 * 1024 * 1024));
	//列式存储，这个大小可以通过控制每次并行的轨迹数目来控制
	CUDA_CALL(cudaMalloc((void**)&this->queryFVGPU, 256 * 1024 * 1024));
	vector<map<int, short>> FVTable;
	FVTable.resize(this->cellNum+1);
	intPair* FVTableOffsetCPU = new intPair[this->cellNum*30000];
	vector<map<int, int>> FVTableOffsetMap;
	FVTableOffsetMap.resize(this->trajNum);

	intPair* FVInfoCPU = new intPair[this->trajNum];
	int cnt = 0;
	for (int i = 1; i <= this->trajNum - 1; i++) {
		map<int, int>::iterator iter;
		FVInfoCPU[i].int_1 = i;
		FVInfoCPU[i].int_2 = cnt;
		for (iter = this->FreqVector[i].begin(); iter != this->FreqVector[i].end(); iter++)
		{
			int cellIdx = iter->first;
			FVTable[cellIdx][i] = iter->second;
			FVTableOffsetMap[i][iter->first] = 0;
		}
	}

	short *FVPtr = (short*)this->FVTableGPU;
	intPair *FVinfoPtr = (intPair*)this->FVinfoGPU;
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	CUDA_CALL(cudaMemcpyAsync(FVinfoPtr, &FVInfoCPU, sizeof(intPair), cudaMemcpyHostToDevice, stream));

	cnt = 0;
	for (int cell = 0; cell <= this->cellNum - 1;cell++)
	{
		map<int, short>::iterator iter;
		for (iter = FVTable[cell].begin(); iter != FVTable[cell].end();iter++)
		{
			CUDA_CALL(cudaMemcpyAsync(FVPtr, &iter->second, sizeof(short), cudaMemcpyHostToDevice, stream));
			//*FVPtr = iter->second;
			FVTableOffsetMap[iter->first][cell] = cnt; //cnt实际上是从FVTable一开始算的offset
			cnt++;
			FVPtr++;
		}
	}
	int FVTableLength = cnt;
	intPair *FVTableOffsetG = (intPair*)this->FVTableOffset;
	for (int i = 1; i <= this->trajNum - 1; i++) {
		map<int, int>::iterator iter;
		for (iter = FVTableOffsetMap[i].begin(); iter != FVTableOffsetMap[i].end();iter++)
		{
			intPair tempPair;
			tempPair.int_1 = iter->first;
			tempPair.int_2 = iter->second;
			CUDA_CALL(cudaMemcpyAsync(FVTableOffsetG, &tempPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
			FVTableOffsetG++;
		}
	}

	//constructing.......

	//int cnt = 0;
	//for (int i = 1; i <= this->trajNum - 1; i++) {
	//	map<int, int>::iterator iter;
	//	intPair tempInfoPair;
	//	tempInfoPair.int_1 = i;
	//	tempInfoPair.int_2 = cnt;
	//	CUDA_CALL(cudaMemcpyAsync(FVinfoPtr, &tempInfoPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
	//	FVinfoPtr++;
	//	for (iter = this->FreqVector[i].begin(); iter != this->FreqVector[i].end();iter++)
	//	{
	//		intPair tempPair;
	//		tempPair.int_1 = iter->first;
	//		tempPair.int_2 = iter->second;
	//		CUDA_CALL(cudaMemcpyAsync(FVPtr, &tempPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
	//		FVPtr++;
	//		cnt++;
	//	}
	//}
	cudaStreamDestroy(stream);
}

int FVTable::formPriorityQueueGPU(priority_queue<FDwithID, vector<FDwithID>, cmp>* queue, map<int, int>* freqVectorQ)
{
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	//传递Query的freqVector到GPU中
	int *queryFV = new int[this->cellNum];
	for (map<int, int>::iterator iter = freqVectorQ->begin(); iter != freqVectorQ->end();iter++)
	{
		queryFV[iter->first] = iter->second;
	}
	int* queryFVGPU = (int*)this->queryFVGPU;
	CUDA_CALL(cudaMemcpyAsync(this->queryFVGPU, queryFV, sizeof(int)*this->cellNum, cudaMemcpyHostToDevice, stream));
	//对于这个查询，调用gpu的kernel执行并行的pruning


	//得到的结果加入到queue中（归并求和）

}

FVTable::FVTable()
{
}


FVTable::~FVTable()
{
}
