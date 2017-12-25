#pragma once
#include "Cell.h"
#include "ConstDefine.h"
#include "MBB.h"
#include "QueryResult.h"
#include "Trajectory.h"
#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include "cudaKernel.h"
#include <map>
#include <bitset>
#include"FVTable.h"

#define MAX_LEVEL (1048576-1)/3

extern Trajectory* tradb;

typedef struct MortonNode {
	int level;
	int nid;
}MortonNode;

class MortonGrid
{
public:
	MBB range;
	float cell_size; //length of a cell
	int cellNum_axis; // 每行/列有几个cell
	int cellnum; //upper(area(grid)/area(cell))，保证能放下所有cell
	Cell* cellPtr; //存储cell的入口
	std::ofstream fout;//文件输出接口
	int totalPointNum; //grid内点个数
	int trajNum;
	int VITURAL_CELL_PARAM;
	std::bitset<size_t(MAX_LEVEL)> isLeaf; // mark whether correspond node is leaf

	SPoint* allPoints;//存储所有点的数组
	Point* allPointsPtrGPU;
	DPoint *allPointsDeltaEncoding;//Delta Encoding后的点
	
	MortonGrid(const MBB& mbb, float val_cell_size, int VITURAL_CELL_PARAM);
	//Range Query on GPU 用
	void *baseAddrRange[2];
	void *stateTableGPU[2];
	RangeQueryStateTable* stateTableRange[2];
	std::map<int, void*> nodeAddrTable[2];
	int stateTableLength[2];
	int nodeAddrTableLength[2];
	int testCnt;
	std::vector<cellBasedTraj> cellBasedTrajectory; //cellbasedtrajectory，二元组：（cell编号数组地址，数组长度）

	//Similarity Query 用
	FVTable freqVectors;

	int addTrajectoryIntoCell(Trajectory &t);
	int WhichCellPointIn(SamplePoint p);
	int addDatasetToGrid(Trajectory* db, int traNum);
	int writeCellsToFile(int* cellNo, int cellNum, std::string file);
	//rangeQuery函数，输入Bounding box，输出轨迹编号和对应顺序下的采样点
	// int rangeQuery(MBB & bound, CPURangeQueryResult * ResultTable, int* resultSetSize);
	//int rangeQueryGPU(MBB & bound, CPURangeQueryResult * ResultTable, int* resultSetSize);
	// int SimilarityQuery(Trajectory &qTra, Trajectory **candTra, int candSize, float *EDRdistance);
	int buildQuadTree(int level, int id);
	MBB generateMBBfromNode(int level, int id);
	//rangeQuery批量
	int rangeQueryBatch(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);
	int rangeQueryBatchMultiThread(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);
	int findMatchNodeInQuadTree(MortonNode node, MBB& bound, std::vector<MortonNode> *cells);
	int rangeQueryBatchGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize, RangeQueryStateTable* stateTableAllocate, int device_idx);
	int rangeQueryBatchMultiGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);
	int findMatchNodeInQuadTreeGPU(MortonNode node, MBB& bound, std::vector<MortonNode> *cells, cudaStream_t stream, int queryID, int device_idx);
	//SimilarityQuery
	int SimilarityQueryBatch(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue);
	int SimilarityQueryBatchCPUParallel(Trajectory *qTra, int queryTrajNum, int *EDRdistance, int kValue);
	int SimilarityMultiThreadHandler(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp>* queryQueue, Trajectory* qTra, int queryTrajNum, std::priority_queue<FDwithID, std::vector<FDwithID>, cmpBig>* EDRCalculated, int kValue, int startQueryIdx);
	int FDCalculateParallelHandeler(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp> *queue, std::map<int, int>* freqVectorQ);
	int SimilarityExecuter(SPoint* queryTra, SPoint** candidateTra, int queryLength, int* candidateLength, int candSize, int *resultArray);
	int SimilarityQueryBatchOnGPU(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue);
	int SimilarityQueryBatchOnMultiGPU(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue);
	MortonGrid();
	~MortonGrid();
};

