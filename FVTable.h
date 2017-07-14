#pragma once
#include "ConstDefine.h"
#include<map>
#include<vector>
#include <queue>

using namespace std;

typedef struct FDwithID {
	int traID;
	int FD;
}FDwithID;

struct cmp {
	bool operator()(FDwithID a, FDwithID b) {
		return(a.FD > b.FD);
	}
};

struct cmpBig {
	bool operator()(FDwithID a, FDwithID b) {
		return(a.FD < b.FD);
	}
};

class FVTable
{
public:
	vector<map<int, int>> FreqVector; //每个map代表一个轨迹
	int trajNum; //总共的轨迹的数量
	int cellNum; //cell的个数
	void *FVTableGPU, *FVinfoGPU,*queryFVGPU,*FVTableOffset;

	int initFVTable(int trajNum, int cellNum);
	int addPointToFVTable(int trajID, int pointNum, int cellID);
	int getCandidate(int lowerBound, int k, map<int, int>* freqVectorQ, int *candidateTrajID, int *candidateNum);
	double calculateFreqDist(int *freqVectorQ, int trajID);
	int findNeighbor(int cellID, int* neighborID);
	int formPriorityQueue(priority_queue<FDwithID, vector<FDwithID>, cmp> *queue, map<int, int>* freqVectorQ);
	// infoFVGPU包含FV在GPU内存储的信息，由于FV在GPU内以一维数组形式存储，每条轨迹的offset需要被保存以便查看。执行完这个函数后，FV table在GPU的地址以及信息保存在成员的变量内。
	int transferFVtoGPU();
	// pruning on GPU
	int formPriorityQueueGPU(priority_queue<FDwithID, vector<FDwithID>, cmp> *queue, map<int, int>* freqVectorQ);

	FVTable();
	~FVTable();
};

