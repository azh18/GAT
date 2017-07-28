#include "STIG.h"
#include "cudaKernel.h"
#include <iostream>
#include <fstream>

using namespace std;
extern void* baseAddrGPU;

STIGNode::STIGNode()
{
}

InternalNode::InternalNode()
{

}

LeafNode::LeafNode()
{

}

STIG::STIG()
{
}

bool cmpBig_SPoint_X(SPoint a, SPoint b) {
	return (a.x > b.x);
}

bool cmpBig_SPoint_Y(SPoint a, SPoint b) {
	return (a.y > b.y);
}

int STIG::initial(int blockSize, int totalDim, Trajectory* db, int trajNum)
{
	this->blockSize = blockSize;
	this->totalDim = totalDim;
	this->maxTid = trajNum;
	int cnt = 0;
	for (int i = 1; i <= trajNum;i++)
	{
		int lengthTraj = db[i].length;
		for (int j = 0; j <= lengthTraj - 1;j++)
		{
			SPoint p;
			p.tID = i;
			p.x = db[i].points[j].lon;
			p.y = db[i].points[j].lat;
			this->allPoints.push_back(p);
			cnt++;
		}
	}
	this->totalPointsNum = cnt;
	if (totalPointsNum <= blockSize)
	{
		this->root = new LeafNode();
		this->createIndex((LeafNode*)this->root, 0, 0, this->totalPointsNum - 1);
	}
	else
	{
		this->root = new InternalNode();
		this->createIndex((InternalNode*)this->root, 0, 0, this->totalPointsNum - 1);
	}
	return 0;
}




int STIG::createIndex(LeafNode* parent, int depth, int startIdx, int endIdx)
{
	parent->leftRange = startIdx;
	parent->rightRange = endIdx;
	parent->depth = depth;
	// generate MBB
	float xmin=999, xmax=-999, ymin=999, ymax=-999;
	for (int i = startIdx; i <= endIdx;i++)
	{
		if (this->allPoints[i].x < xmin)
			xmin = this->allPoints[i].x;
		if (this->allPoints[i].x > xmax)
			xmax = this->allPoints[i].x;
		if (this->allPoints[i].y < ymin)
			ymin = this->allPoints[i].y;
		if (this->allPoints[i].y > ymax)
			ymax = this->allPoints[i].y;
	}
	MBB mbb(xmin, ymin, xmax, ymax);
	parent->boundingBox = mbb;
	//printf("[%d,%d]:%f,%f,%f,%f------", parent->leftRange, parent->rightRange, mbb.xmin, mbb.xmax, mbb.ymin, mbb.ymax);
	return 0;
}


int STIG::createIndex(InternalNode* parent, int depth, int startIdx, int endIdx)
{
	parent->depth = depth;
	if (depth % 2 == 0)
	{
		std::sort(&this->allPoints[startIdx], &this->allPoints[endIdx], cmpBig_SPoint_X);
		int medianIdx = (startIdx + endIdx) / 2;
		float medianX = this->allPoints[medianIdx].x;
		parent->dim = depth % 2;
		parent->medium = medianX;
		// left child
		if (medianIdx - startIdx <= this->blockSize)
		{
			parent->left = new LeafNode();
			this->createIndex((LeafNode*)parent->left, depth + 1, startIdx, medianIdx);
			parent->leftIsLeaf = TRUE;
		}
		else
		{
			parent->left = new InternalNode();
			this->createIndex((InternalNode*)parent->left, depth + 1, startIdx, medianIdx-1);
			parent->leftIsLeaf = FALSE;
		}
		//right child
		if(endIdx - medianIdx + 1<=this->blockSize)
		{
			parent->right = new LeafNode();
			this->createIndex((LeafNode*)parent->right, depth + 1, medianIdx, endIdx);
			parent->rightIsLeaf = TRUE;
		}
		else
		{
			parent->right = new InternalNode();
			this->createIndex((InternalNode*)parent->right, depth + 1, medianIdx, endIdx);
			parent->rightIsLeaf = FALSE;
		}
	}
	else
	{
		std::sort(&this->allPoints[startIdx], &this->allPoints[endIdx], cmpBig_SPoint_Y);
		int medianIdx = (startIdx + endIdx) / 2;
		float medianY = this->allPoints[medianIdx].y;
		parent->dim = depth % 2;
		parent->medium = medianY;
		// left child
		if (medianIdx - startIdx <= this->blockSize)
		{
			parent->left = new LeafNode();
			this->createIndex((LeafNode*)parent->left, depth + 1, startIdx, medianIdx);
			parent->leftIsLeaf = TRUE;
		}
		else
		{
			parent->left = new InternalNode();
			this->createIndex((InternalNode*)parent->left, depth + 1, startIdx, medianIdx - 1);
			parent->leftIsLeaf = FALSE;
		}
		//right child
		if (endIdx - medianIdx + 1 <= this->blockSize)
		{
			parent->right = new LeafNode();
			this->createIndex((LeafNode*)parent->right, depth + 1, medianIdx, endIdx);
			parent->rightIsLeaf = TRUE;
		}
		else
		{
			parent->right = new InternalNode();
			this->createIndex((InternalNode*)parent->right, depth + 1, medianIdx, endIdx);
			parent->rightIsLeaf = FALSE;
		}
	}
	return 0;
}



int STIG::destroyIndex()
{
	
}



int STIG::searchTree(MBB queryMBB, std::vector<STIGBlock>* allCandBlocks)
{
	if(this->totalPointsNum<=this->blockSize)
	{
		//root is leaf
		LeafNode* root = (LeafNode*)this->root;
		if(root->boundingBox.intersect(queryMBB))
		{
			STIGBlock b;
			b.startIdx = root->leftRange;
			b.endIdx = root->rightRange;
			allCandBlocks->push_back(b);
		}
		return 0;
	}
	else
	{
		this->searchNode(queryMBB, allCandBlocks, (InternalNode*)root);
	}
	return 0;
}

int STIG::searchNode(MBB queryMBB, std::vector<STIGBlock>* allCandBlocks, InternalNode* node)
{
	int dim = node->depth % 2;
	if(dim == 0)
	{
		if(intersectBlock(queryMBB.xmin,queryMBB.xmax,0,node->medium))
		{
			if (node->leftIsLeaf)
				searchNode(queryMBB, allCandBlocks, (LeafNode*)node->left);
			else
				searchNode(queryMBB, allCandBlocks, (InternalNode*)node->left);
		}
		if(intersectBlock(queryMBB.xmin,queryMBB.xmax,node->medium,180))
		{
			if (node->rightIsLeaf)
				searchNode(queryMBB, allCandBlocks, (LeafNode*)node->right);
			else
				searchNode(queryMBB, allCandBlocks, (InternalNode*)node->right);
		}
	}
	else
	{
		if (intersectBlock(queryMBB.ymin, queryMBB.ymax, 0, node->medium))
		{
			if (node->leftIsLeaf)
				searchNode(queryMBB, allCandBlocks, (LeafNode*)node->left);
			else
				searchNode(queryMBB, allCandBlocks, (InternalNode*)node->left);
		}
		if (intersectBlock(queryMBB.xmin, queryMBB.xmax, node->medium, 180))
		{
			if (node->rightIsLeaf)
				searchNode(queryMBB, allCandBlocks, (LeafNode*)node->right);
			else
				searchNode(queryMBB, allCandBlocks, (InternalNode*)node->right);
		}
	}
	return 0;
}

int STIG::searchNode(MBB queryMBB, std::vector<STIGBlock>* allCandBlocks, LeafNode* node)
{
	if (node->boundingBox.intersect(queryMBB))
	{
		STIGBlock b;
		b.startIdx = node->leftRange;
		b.endIdx = node->rightRange;
		allCandBlocks->push_back(b);
		//printf("Result: [%d,%d]:%f,%f,%f,%f\n", node->leftRange, node->rightRange, node->boundingBox.xmin, node->boundingBox.xmax, node->boundingBox.ymin, node->boundingBox.ymax);
	}
	return 0;
}

int STIG::rangeQueryGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize)
{
	// search all blocks should be sent to GPU
	cudaStream_t stream;
	cudaStreamCreate(&stream);

	void* baseAddr = baseAddrGPU;
	void* allPointsInGPU = baseAddr;
	vector<RangeQueryStateTable> stateTable;
	int cntTask = 0;
	for (int i = 0; i <= rangeNum - 1; i++) {
		vector<STIGBlock> allCandBlocks;
		this->searchTree(bounds[i], &allCandBlocks);
		for (int j = 0; j <= allCandBlocks.size() - 1;j++)
		{
			RangeQueryStateTable stateTableCPU;
			int startIdx = allCandBlocks[j].startIdx;
			int endIdx = allCandBlocks[j].endIdx;
			int candPointsNum = endIdx - startIdx + 1;
			//SPoint pp[1000];
			//memcpy(pp, &this->allPoints[startIdx], sizeof(SPoint)*candPointsNum);
			//printf("startIdx: %d\n", startIdx);
			CUDA_CALL(cudaMemcpyAsync(baseAddr, (void*)(&this->allPoints[startIdx]), sizeof(SPoint)*candPointsNum, cudaMemcpyHostToDevice, stream));
			stateTableCPU.ptr = baseAddr;
			stateTableCPU.xmin = bounds[i].xmin;
			stateTableCPU.xmax = bounds[i].xmax;
			stateTableCPU.ymin = bounds[i].ymin;
			stateTableCPU.ymax = bounds[i].ymax;
			stateTableCPU.queryID = i;
			stateTableCPU.candidatePointNum = candPointsNum;
			stateTableCPU.startIdxInAllPoints = 0;
			baseAddr = (void*)((SPoint*)baseAddr + candPointsNum);
			stateTable.push_back(stateTableCPU);
			cntTask++;
		}
	}
	// allocate GPU, using stateTable
	RangeQueryStateTable* stateTableGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)&stateTableGPU, sizeof(RangeQueryStateTable)*stateTable.size()));
	CUDA_CALL(cudaMemcpyAsync(stateTableGPU, (void*)(&stateTable[0]), sizeof(RangeQueryStateTable)*stateTable.size(), cudaMemcpyHostToDevice, stream));
	//RangeQueryStateTable* stateTableGPU = (RangeQueryStateTable*)baseAddr;
	uint8_t *resultsReturned = new uint8_t[rangeNum*(this->maxTid+1)];
	memset(resultsReturned, 0, rangeNum*(this->maxTid+1));
	cudaRangeQuerySTIGHandler(stateTableGPU, stateTable.size(), resultsReturned, (this->maxTid+1), rangeNum, stream);
	//shanhou
	ofstream fp("queryResult(STIG).txt",ios_base::out);
	for (int jobID = 0; jobID <= rangeNum - 1; jobID++) {
		for (int traID = 1; traID <= (this->maxTid); traID++) {
			if (resultsReturned[jobID*(this->maxTid+1) + traID] == 1) {
				fp << "job " << jobID << "find" << traID << endl;
			}
		}
	}
	cudaStreamDestroy(stream);
	return 0;
}

bool STIG::intersectBlock(int amin,int amax,int bmin,int bmax)
{
	if (amax < bmin)
		return FALSE;
	if (amin > bmax)
		return FALSE;
	return TRUE;
}
