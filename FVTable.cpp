#include "FVTable.h"
#include <queue>
#include "cudaKernel.h"

using namespace std;

int FVTable::initFVTable(int trajNum,int cellNum)
{
	this->cellNum = cellNum;
	this->trajNum = trajNum;
	this->FreqVector.resize(trajNum+1);
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
	for (int i = 1; i <= this->trajNum; i++) {
		//对于该轨迹，计算与查询的FV的FD
		//首先，计算两个vector的对应元素减法
		//原始轨迹减去查询轨迹
		//time1.start();
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
		//time1.stop();
		//printf("prun time 1:%f\n", time1.elapse());
		//time1.start();
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
		//time1.stop();
		//printf("prun time 2:%f\n", time1.elapse());
		//time1.start();
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
		//time1.stop();
		//printf("prun time 3:%f\n", time1.elapse());
		//time1.start();
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
		//time1.stop();
		//printf("prun time 4:%f\n", time1.elapse());
		//time1.start();
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
		//time1.stop();
		//printf("prun time 5:%f\n", time1.elapse());
	}
	return 0;
}

int FVTable::transferFVtoGPU()
{
	cudaStream_t stream;
	cudaStreamCreate(&stream);
#ifdef NOT_COLUMN_ORIENTED
	CUDA_CALL(cudaMalloc((void**)&this->FVinfoGPU, 16 * 1024 * 1024));
	CUDA_CALL(cudaMalloc((void**)&this->FVTableOffset, 256 * 1024 * 1024));
	CUDA_CALL(cudaMalloc((void**)&this->FVTableGPU, 256 * 1024 * 1024));
	CUDA_CALL(cudaMalloc((void**)&this->queryFVGPU, this->cellNum*sizeof(short)*N_BATCH_QUERY));
	CUDA_CALL(cudaMalloc((void**)&this->FDresultsGPU, N_BATCH_QUERY * sizeof(short)));
	intPair* FVinfoPtr = (intPair*)this->FVinfoGPU;
	intPair* FVPtr = (intPair*)this->FVTableGPU;
	int cnt = 0;
	for (int i = 1; i <= this->trajNum; i++) {
		map<int, int>::iterator iter;
		intPair tempInfoPair;
		tempInfoPair.int_1 = i;
		tempInfoPair.int_2 = cnt;
		CUDA_CALL(cudaMemcpyAsync(FVinfoPtr, &tempInfoPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
		FVinfoPtr++;
		for (iter = this->FreqVector[i].begin(); iter != this->FreqVector[i].end(); iter++)
		{
			intPair tempPair;
			tempPair.int_1 = iter->first;
			tempPair.int_2 = iter->second;
			CUDA_CALL(cudaMemcpyAsync(FVPtr, &tempPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
			FVPtr++;
			cnt++;
		}
	}
	this->nonZeroFVNum = cnt;

#else

	
	CUDA_CALL(cudaMalloc((void**)&this->FVinfoGPU, sizeof(intPair) * 80000));// 保存trajID、length、在FVTable中的offset
	//CUDA_CALL(cudaMalloc((void**)&this->FVTableOffset, sizeof(intPair) * 80000 * 1024)); // 保存cellID、在FVTable中的offset
	CUDA_CALL(cudaMalloc((void**)&this->FVTableGPU, sizeof(intPair) * 80000 * 1024)); // 保存在数据库中的（cellID,freq）数据
	CUDA_CALL(cudaMalloc((void**)&this->FDresultsGPU, N_BATCH_QUERY * sizeof(int)));
	intPair* FVinfoPtr = (intPair*)this->FVinfoGPU;
	intPair* FVPtr = (intPair*)this->FVTableGPU;
	//cnt记录在所有轨迹中，某个轨迹的编号范围
	int cnt = 0;
	int maxTrajCellLength = 0;
	for (int i = 1; i <= this->trajNum; i++) {
		map<int, int>::iterator iter;
		intPair tempInfoPair;
		tempInfoPair.int_1 = i;
		tempInfoPair.int_2 = cnt;
		CUDA_CALL(cudaMemcpyAsync(FVinfoPtr, &tempInfoPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
		FVinfoPtr++;
		for (iter = this->FreqVector[i].begin(); iter != this->FreqVector[i].end(); iter++)
		{
			intPair tempPair;
			tempPair.int_1 = iter->first;
			tempPair.int_2 = iter->second;
			CUDA_CALL(cudaMemcpyAsync(FVPtr, &tempPair, sizeof(intPair), cudaMemcpyHostToDevice, stream));
			FVPtr++;
			cnt++;
		}
		if (this->FreqVector[i].size() > maxTrajCellLength)
			maxTrajCellLength = this->FreqVector[i].size();
	}
	this->SubbedArrayJump = 2*maxTrajCellLength;
	CUDA_CALL(cudaMalloc((void**)&this->queryFVGPU, (maxTrajCellLength)*sizeof(intPair)));// 存储query轨迹的FV
	CUDA_CALL(cudaMalloc((void**)&this->SubbedArrayGPU, sizeof(intPair) * (this->SubbedArrayJump) * N_BATCH_QUERY));// 预留给GPU用的用于计算的空间
	//printf("%d", maxTrajCellLength);
	CUDA_CALL(cudaMalloc((void**)&this->SubbedArrayOffsetGPU, sizeof(intPair) * N_BATCH_QUERY));// 预留给GPU用的用于计算的空间
	this->nonZeroFVNum = cnt;

#endif

	cudaStreamDestroy(stream);
	return 0;
}

int FVTable::formPriorityQueueGPU(priority_queue<FDwithID, vector<FDwithID>, cmp>* queue, map<int, int>* freqVectorQ)
{
	cudaStream_t stream;
	cudaStreamCreate(&stream);
#ifdef NOT_COLUMN_ORIENTED
	short *queryFVCPU = (short*)malloc(sizeof(short)*this->cellNum*N_BATCH_QUERY);
	for (map<int, int>::iterator iter = freqVectorQ->begin(); iter != freqVectorQ->end(); iter++)
	{
		for (int line = 0; line <= N_BATCH_QUERY - 1; line++)
			queryFVCPU[line*this->cellNum + iter->first] = iter->second;
	}
	short* queryFVGPU = (short*)this->queryFVGPU;

	int candidateTotalNum = this->trajNum;
	for (int i = 1; i <= trajNum; i += N_BATCH_QUERY)
	{
		//一次计算taskNum个FD，由于GPU内存限制
		int taskNum = N_BATCH_QUERY;
		if (i + N_BATCH_QUERY > trajNum)
			taskNum = trajNum - i + 1;
		CUDA_CALL(cudaMemcpyAsync(queryFVGPU, queryFVCPU, sizeof(short)*this->cellNum*taskNum, cudaMemcpyHostToDevice, stream));
		//对于这个查询，调用gpu的kernel执行并行的pruning (注意传入pitch)
		//这里的trajIdx从0开始,checkNum是指待检查的轨迹的条数（block数量）
		Similarity_Pruning_Handler((short*)this->queryFVGPU, (intPair*)this->FVinfoGPU, (intPair*)this->FVTableGPU, i - 1, taskNum, this->cellNum, this->trajNum, this->nonZeroFVNum, (short*)this->FDresultsGPU, stream);
		short* resultsTemp = new short[taskNum];
		CUDA_CALL(cudaMemcpyAsync(resultsTemp, (short*)this->FDresultsGPU, sizeof(short)*taskNum, cudaMemcpyDeviceToHost, stream));
		//得到的结果加入到queue中（归并求和）
		for (int item = i; item < i + taskNum; item++) {
			FDwithID fd;
			fd.traID = item;
			fd.FD = resultsTemp[item-i-1];
			queue->push(fd);
		}
		delete[] resultsTemp;
	}
	free(queryFVCPU);
#else
	intPair *queryFVCPU = (intPair*)malloc(freqVectorQ->size());
	int queryCellLength = 0;
	for (map<int, int>::iterator iter = freqVectorQ->begin(); iter != freqVectorQ->end(); iter++)
	{
		queryFVCPU[queryCellLength].int_1 = iter->first;
		queryFVCPU[queryCellLength].int_2 = iter->second;
		queryCellLength++;
	}
	intPair* queryFVGPU = (intPair*)this->queryFVGPU;

	int candidateTotalNum = this->trajNum;
	for (int i = 1; i <= trajNum; i += N_BATCH_QUERY)
	{
		//一次计算taskNum个FD，由于GPU内存限制
		int taskNum = N_BATCH_QUERY;
		if (i + N_BATCH_QUERY > trajNum)
			taskNum = trajNum - i + 1;
		CUDA_CALL(cudaMemcpyAsync(queryFVGPU, queryFVCPU, sizeof(intPair)*queryCellLength, cudaMemcpyHostToDevice, stream));
		//对于这个查询，调用gpu的kernel执行并行的pruning (注意传入pitch)
		//这里的trajIdx从0开始,checkNum是指待检查的轨迹的条数（block数量）
		Similarity_Pruning_Handler((intPair*)this->queryFVGPU, (intPair*)this->FVinfoGPU, (intPair*)this->FVTableGPU,(intPair*)this->SubbedArrayGPU, (intPair*)SubbedArrayOffsetGPU,this->SubbedArrayJump, queryCellLength, i - 1, taskNum, this->cellNum, this->trajNum, this->nonZeroFVNum, (short*)this->FDresultsGPU, stream);
		short* resultsTemp = new short[taskNum];
		CUDA_CALL(cudaMemcpyAsync(resultsTemp, (short*)this->FDresultsGPU, sizeof(short)*taskNum, cudaMemcpyDeviceToHost, stream));
		//得到的结果加入到queue中（归并求和）
		for (int item = i; item < i + taskNum; item++) {
			FDwithID fd;
			fd.traID = item;
			fd.FD = resultsTemp[item - i - 1];
			queue->push(fd);
		}
		delete[] resultsTemp;
	}
	free(queryFVCPU);

#endif
	
	return 0;
}

FVTable::FVTable()
{
}


FVTable::~FVTable()
{
}
