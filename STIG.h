#pragma once

#include "ConstDefine.h"
#include "MBB.h"
#include "Trajectory.h"
#include <algorithm>
#include <vector>

typedef struct STIGBlock{
	int startIdx;
	int endIdx;
}STIGBlock;

class STIGNode
{
public:
	STIGNode();
	int depth;
};

class InternalNode: public STIGNode
{
public:
	InternalNode();
	int dim;
	float medium;
	bool leftIsLeaf;
	bool rightIsLeaf;
	void* left;
	void* right;
};

class LeafNode: public STIGNode
{
public:
	LeafNode();
	MBB boundingBox;
	int leftRange, rightRange;
};

class STIG
{
public:
	STIG();
	// 多GPU查询用
	void* baseAddrGPU[2];
	void* stateTableGPU[2];
	// 多GPU查询用

	void *root;
	int totalDim;
	int totalPointsNum;
	std::vector<SPoint> allPoints;
	int depth;
	int blockSize;
	int maxTid;
	int initial(int blockSize, int totalDim, Trajectory* db, int trajNum);
	int createIndex(LeafNode* parent, int depth, int startIdx, int endIdx);
	int createIndex(InternalNode* parent, int depth, int startIdx, int endIdx);
	int destroyIndex();
	int searchTree(MBB queryMBB, std::vector<STIGBlock> *allCandBlocks);
	int searchNode(MBB queryMBB, std::vector<STIGBlock> *allCandBlocks, InternalNode* node);
	int searchNode(MBB queryMBB, std::vector<STIGBlock> *allCandBlocks, LeafNode* node);
	int rangeQueryGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize, int device_idx);
	static bool intersectBlock(float amin, float amax, float bmin, float bmax);
};