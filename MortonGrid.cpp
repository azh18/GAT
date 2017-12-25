#include "MortonGrid.h"

using namespace std;

extern int getIdxFromXY(int x, int y);


inline int findIdxInBitmap(int level, int id) {
	if (level <= 0)
		return 0;
	else
		return (int(pow(4, level) + 0.5) - 1) / 3 + id;
}
MortonGrid::MortonGrid()
{
}


MortonGrid::~MortonGrid()
{
}

MortonGrid::MortonGrid(const MBB & mbb, float val_cell_size, int VITURAL_CELL_PARAM)
{
	this->VITURAL_CELL_PARAM = VITURAL_CELL_PARAM;
	range = mbb;
	cell_size = val_cell_size;
	//貌似只需要用一个维度就行了，因为规定好了必须是2*2,4*4，……
	int divideNumOnX = (int)((mbb.xmax - mbb.xmin) / val_cell_size) + 1; //最少要用多少个cell
	int divideNumOnY = (int)((mbb.ymax - mbb.ymin) / val_cell_size) + 1;
	int maxValue = max(divideNumOnX, divideNumOnY);
	//找到最佳的长宽
	cellNum_axis = maxValue >> (int(log2(maxValue))) << (int(log2(maxValue)) + 1);
	cellnum = cellNum_axis * cellNum_axis;
	cellPtr = new Cell[cellnum];
	this->testCnt = 0;
	//由于满足正方形需要，向xmax、ymin方向扩展range

	//注意cell编号是从(xmin,ymax)开始的，而不是(xmin,ymin)
	//Z字形编码
	for (int i = 0; i <= cellNum_axis - 1; i++)
	{
		for (int j = 0; j <= cellNum_axis - 1; j++)
		{
			int cell_idx = getIdxFromXY(j, i);
			cellPtr[cell_idx].initial(i, j, MBB(range.xmin + cell_size * j, range.ymax - cell_size * (i + 1), range.xmin + cell_size * (j + 1), range.ymax - cell_size * (i)));
		}
	}
}

int MortonGrid::addTrajectoryIntoCell(Trajectory & t)
{
	if (t.length == 0)
		return 1;//空轨迹
	SamplePoint p = t.points[0];
	int lastCellNo = WhichCellPointIn(p);
	int lastCellStartIdx = 0;
	int nowCellNo;
	//cell based traj生成，记得转换后free！
	vector<int>* tempCellBasedTraj = new vector<int>;
	tempCellBasedTraj->reserve(1048577);
	int tempCellNum = 0;
	for (int i = 0; i <= t.length - 1; i++)
	{
		p = t.points[i];
		nowCellNo = WhichCellPointIn(p);
		if (i == t.length - 1)
		{
			//到最后一条，发现这个cell也是上个cell就是最后一个cell了，添加之
			if (lastCellNo == nowCellNo)
			{
				tempCellNum++;
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[nowCellNo].addSubTra(t.tid, lastCellStartIdx, i, i - lastCellStartIdx + 1);
				int vituralCellNo = nowCellNo >> VITURAL_CELL_PARAM; //虚格子
				this->freqVectors.addPointToFVTable(t.tid, i - lastCellStartIdx + 1, vituralCellNo);
			}
			//否则，上一个和这个cell都要添加
			else
			{
				tempCellNum += 2;
				tempCellBasedTraj->push_back(lastCellNo);
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				int vituralCellNo = lastCellNo >> VITURAL_CELL_PARAM; //虚格子
				this->freqVectors.addPointToFVTable(t.tid, i - 1 - lastCellStartIdx + 1, vituralCellNo);
				cellPtr[nowCellNo].addSubTra(t.tid, i, i, 1);
				vituralCellNo = nowCellNo >> VITURAL_CELL_PARAM;
				this->freqVectors.addPointToFVTable(t.tid, 1, vituralCellNo);
			}
		}
		else
		{
			if (lastCellNo == nowCellNo)
				continue;
			else
			{
				// 终结一条子轨迹，开始下一条子轨迹
				//cellTra里面加一条
				tempCellNum++;
				tempCellBasedTraj->push_back(lastCellNo);
				//SubTra添加
				//printf("cell:%d\n", lastCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				int vituralCellNo = lastCellNo >> VITURAL_CELL_PARAM; //虚格子
				this->freqVectors.addPointToFVTable(t.tid, i - 1 - lastCellStartIdx + 1, vituralCellNo);
				lastCellNo = nowCellNo;
				lastCellStartIdx = i;
			}
		}
	}
	this->cellBasedTrajectory[t.tid].length = tempCellNum;
	this->cellBasedTrajectory[t.tid].cellNo = (int*)malloc(sizeof(int) * tempCellNum);
	if (this->cellBasedTrajectory[t.tid].cellNo == NULL) throw("alloc error");
	for (int i = 0; i <= tempCellNum - 1; i++)
	{
		this->cellBasedTrajectory[t.tid].cellNo[i] = tempCellBasedTraj->at(i);
	}
	this->cellBasedTrajectory[t.tid].trajLength = t.length;
	delete tempCellBasedTraj;
	return 0;
	return 0;
}

int MortonGrid::WhichCellPointIn(SamplePoint p)
{
	//注意cell编号是从(xmin,ymax)开始的，而不是(xmin,ymin)
	int row = (int)((range.ymax - p.lat) / cell_size); //从0开始
	int col = (int)((p.lon - range.xmin) / cell_size); //从0开始
	return getIdxFromXY(col, row);
	return 0;
}

int MortonGrid::addDatasetToGrid(Trajectory * db, int traNum)
{
	this->trajNum = traNum;
	//生成frequency vector
	this->freqVectors.initFVTable(traNum, (this->cellnum) >> 2); //注意这里是虚格子的个数
																 //注意，轨迹编号从1开始
	this->cellBasedTrajectory.resize(traNum + 1); //扩大cellbasedtraj的规模，加轨迹的时候可以直接用
	int pointCount = 0;
	for (int i = 1; i <= traNum; i++)
	{
		addTrajectoryIntoCell(db[i]);
	}
	for (int i = 0; i <= cellnum - 1; i++)
	{
		cellPtr[i].buildSubTraTable();
		pointCount += cellPtr[i].totalPointNum;
	}
	this->totalPointNum = pointCount;
	//链表变成了数组
	//subTraTable仅仅是记录了子轨迹（起始offset、终止offset、Tid）

	//建立Quadtree，自顶向下建立，分割节点使所有节点包含点的个数小于MAXPOINTINNODE
	// this->root = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));
	this->buildQuadTree(0, 0);

	//转化为cell连续存储
	//此处连续存储是指同一cell内的采样点存储在一起，有利于rangeQuery，但不利于similarity query
	//similarity组装轨迹的时候，可以先记录当前是第几个subtra，找轨迹的时候从这个往后找，避免tid重复存在的问题
	this->allPoints = (SPoint*)malloc(sizeof(SPoint) * (this->totalPointNum));
	pointCount = 0;


	for (int i = 0; i <= cellnum - 1; i++)
	{
		cellPtr[i].pointRangeStart = pointCount;
		for (int j = 0; j <= cellPtr[i].subTraNum - 1; j++)
		{
			//for each subTra, add Points to AllPoints
			cellPtr[i].subTraTable[j].idxInAllPointsArray = pointCount;
			for (int k = cellPtr[i].subTraTable[j].startpID; k <= cellPtr[i].subTraTable[j].endpID; k++)
			{
				allPoints[pointCount].tID = cellPtr[i].subTraTable[j].traID;
				allPoints[pointCount].x = tradb[allPoints[pointCount].tID].points[k].lon;
				allPoints[pointCount].y = tradb[allPoints[pointCount].tID].points[k].lat;
				//allPoints[pointCount].time = tradb[allPoints[pointCount].tID].points[k].time;
				pointCount++;
			}
		}
		cellPtr[i].pointRangeEnd = pointCount - 1;
		if (cellPtr[i].pointRangeEnd - cellPtr[i].pointRangeStart + 1 != cellPtr[i].totalPointNum)
			cerr << "Grid.cpp: something wrong in total point statistic" << endl;
	}
	//完善CellBasedTrajectory，主要是添加idx信息
	//对每条轨迹分配cellBasedTrajectory后两个信息的内存：起始位置和采样点数目
	for (int i = 1; i <= this->trajNum; i++)
	{
		this->cellBasedTrajectory[i].startIdx = (int*)malloc(sizeof(int) * this->cellBasedTrajectory[i].length);
		this->cellBasedTrajectory[i].numOfPointInCell = (short*)malloc(sizeof(short) * this->cellBasedTrajectory[i].length);
		int* tempCntForTraj = (int*)malloc(sizeof(int) * this->cellnum);
		memset(tempCntForTraj, 0, sizeof(int) * this->cellnum);
		for (int cellidx = 0; cellidx <= this->cellBasedTrajectory[i].length - 1; cellidx++)
		{
			int nowCellID = this->cellBasedTrajectory[i].cellNo[cellidx];
			int j, cnt;
			for (j = 0, cnt = 0; cnt <= tempCntForTraj[nowCellID]; j++)
			{
				if (this->cellPtr[nowCellID].subTraTable[j].traID == i)
				{
					cnt++;
				}
			}
			j--;
			//choose j; //选择第j项subTra，但是这里面还是没有起始点的信息，考虑回来把point的信息加进subtraTable。。。反正最后要把subTraTable给free掉
			this->cellBasedTrajectory[i].startIdx[cellidx] = this->cellPtr[nowCellID].subTraTable[j].idxInAllPointsArray;
			this->cellBasedTrajectory[i].numOfPointInCell[cellidx] = this->cellPtr[nowCellID].subTraTable[j].numOfPoint;
			tempCntForTraj[nowCellID]++;
		}
		free(tempCntForTraj);
	}
	return 0;
}

int MortonGrid::writeCellsToFile(int * cellNo, int cellNum, std::string file)
{
	fout.open(file, ios_base::out);
	for (int i = 0; i <= cellNum - 1; i++)
	{
		int outCellIdx = cellNo[i];
		cout << outCellIdx << ": " << "[" << cellPtr[outCellIdx].mbb.xmin << "," << cellPtr[outCellIdx].mbb.xmax << "," << cellPtr[outCellIdx].mbb.ymin << "," << cellPtr[outCellIdx].mbb.ymax << "]" << endl;
		for (int j = 0; j <= cellPtr[outCellIdx].subTraNum - 1; j++)
		{
			int tid = cellPtr[outCellIdx].subTraTable[j].traID;
			int startpid = cellPtr[outCellIdx].subTraTable[j].startpID;
			int endpid = cellPtr[outCellIdx].subTraTable[j].endpID;
			for (int k = startpid; k <= endpid; k++)
			{
				cout << tradb[tid].points[k].lat << "," << tradb[tid].points[k].lon << ";";
			}
			cout << endl;
		}
	}
	return 0;
}

int MortonGrid::buildQuadTree(int level, int id)
{
	int totalLevel = int(log2(this->cellnum) / log2(4));
	int totalPoints = 0;
	for (int i = id * int(pow(4, (totalLevel - level))); i <= (id + 1) * int(pow(4, (totalLevel - level))) - 1; i++)
	{
		totalPoints += this->cellPtr[i].totalPointNum;
	}
	if ((totalPoints < MAXPOINTINNODE) || (level == totalLevel))
	{
		int idxInBitmap = findIdxInBitmap(level, id);
		this->isLeaf[idxInBitmap] = 1;
		return 0;
	}
	else
	{
		this->buildQuadTree(level + 1, id << 2);
		this->buildQuadTree(level + 1, (id << 2) + 1);
		this->buildQuadTree(level + 1, (id << 2) + 2);
		this->buildQuadTree(level + 1, (id << 2) + 3);
		return 0;
	}
	return 0;
}

MBB MortonGrid::generateMBBfromNode(int level, int id)
{
	int totalLevel = int(log2(this->cellnum) / log2(4));
	return MBB(this->cellPtr[id * int(pow(4, (totalLevel - level)))].mbb.xmin, this->cellPtr[(id + 1) * int(pow(4, (totalLevel - level))) - 1].mbb.ymin, this->cellPtr[(id + 1) * int(pow(4, (totalLevel - level))) - 1].mbb.xmax, this->cellPtr[id * int(pow(4, (totalLevel - level)))].mbb.ymax);
}


int MortonGrid::rangeQueryBatch(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize)
{
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		ResultTable[i].resize(this->trajNum + 1);
	}
#ifdef CHECK_CORRECT
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		for (int j = 0; j <= this->trajNum + 1; j++)
		{
			ResultTable[i][j] = 0;
		}
	}
#endif

	int totalLevel = int(log2(this->cellnum) / log2(4));
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		//int candidateNodeNum = 0;
		vector<MortonNode> cells;
		MortonNode root;
		root.level = 0;
		root.nid = 0;
		findMatchNodeInQuadTree(root, bounds[i], &cells);
		//printf("%d", cells.size());
		for (vector<MortonNode>::iterator iterV = cells.begin(); iterV != cells.end(); iterV++)
		{
			int nodeID = iterV->nid;
			int nodeLevel = iterV->level;
			int firstCellID = nodeID * int(pow(4, (totalLevel - nodeLevel))+0.5);
			int lastCellID = (nodeID + 1) * int(pow(4, (totalLevel - nodeLevel))+0.5) - 1;
			for (int cellID = firstCellID; cellID <= lastCellID; cellID++)
			{
				int anchorX = this->cellPtr[cellID].anchorPointX;
				int anchorY = this->cellPtr[cellID].anchorPointY;
				for (int idx = this->cellPtr[cellID].pointRangeStart; idx <= this->cellPtr[cellID].pointRangeEnd; idx++)
				{
					//compress
					//float realX = float(allPointsDeltaEncoding[idx].x + anchorX) / 1000000;
					//float realY = float(allPointsDeltaEncoding[idx].y + anchorY) / 1000000;
					// no compress
					float realX = allPoints[idx].x;
					float realY = allPoints[idx].y;
					int tID = allPoints[idx].tID;
					if (bounds[i].pInBox(realX, realY))
					{
						ResultTable[i][tID] = TRUE;
					}
				}
			}
		}
	}

	//for (int jobID = 0; jobID <= rangeNum - 1; jobID++)
	//{
	//	for (int traID = 1; traID <= (this->trajNum); traID++)
	//	{
	//		if (resultsReturned[jobID * (this->trajNum + 1) + traID] == 1)
	//		{
	//			out << "job " << jobID << "find" << traID << endl;
	//		}
	//	}
	//}
	return 0;
}

int MortonGrid::rangeQueryBatchMultiThread(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize)
{
	vector<thread> threads_RQ;
	const int device_num = 16;
	int nowidx = 0;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		int boundNum = rangeNum / device_num + 1;
		if (nowidx + boundNum >= rangeNum)
			boundNum = rangeNum - nowidx;
		threads_RQ.push_back(thread(std::mem_fn(&MortonGrid::rangeQueryBatch), this, &bounds[nowidx], boundNum, &ResultTable[nowidx], resultSetSize));
		nowidx += boundNum;
	}
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));

	return 0;
}

int MortonGrid::findMatchNodeInQuadTree(MortonNode node, MBB & bound, std::vector<MortonNode>* cells)
{
	if (this->isLeaf[findIdxInBitmap(node.level,node.nid)])
	{
		cells->push_back(node);
	}
	else
	{
		MBB UL = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 0);
		if (bound.intersect(UL)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 0;
			findMatchNodeInQuadTree(tempNode, bound, cells);
		}
		MBB UR = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 1);
		if (bound.intersect(UR)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 1;
			findMatchNodeInQuadTree(tempNode, bound, cells);
		}
		MBB DL = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 2);
		if (bound.intersect(DL)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 2;
			findMatchNodeInQuadTree(tempNode, bound, cells);
		}
		MBB DR = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 3);
		if (bound.intersect(DR)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 3;
			findMatchNodeInQuadTree(tempNode, bound, cells);
		}
	}
	return 0;
}

int MortonGrid::rangeQueryBatchGPU(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize, RangeQueryStateTable * stateTableAllocate, int device_idx)
{
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		ResultTable[i].resize(this->trajNum + 1);
	}
#ifdef CHECK_CORRECT
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		for (int j = 0; j <= this->trajNum; j++)
		{
			ResultTable[i][j] = 0;
		}
	}
#endif
	// 分配GPU内存
	//MyTimer timer;
	// 参数随便设置的，可以再调
	//timer.start();
	CUDA_CALL(cudaSetDevice(device_idx));
	this->stateTableRange[device_idx] = stateTableAllocate;
	this->stateTableLength[device_idx] = 0;
	this->nodeAddrTableLength[device_idx] = 0;
	// for each query, generate the nodes:
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		MortonNode root;
		root.level = 0;
		root.nid = 0;
		findMatchNodeInQuadTreeGPU(root, bounds[i], NULL, stream, i, device_idx);
	}
	//printf("StateTableLength:%d",this->stateTableLength);
	//stateTable中点的数目的最大值
	FILE* testFile = fopen("freq.txt", "w+");
	int maxPointNum = 0;
	for (int i = 0; i <= this->stateTableLength[device_idx] - 1; i++)
	{
		if (stateTableAllocate[i].candidatePointNum > maxPointNum)
			maxPointNum = stateTableAllocate[i].candidatePointNum;
		//fprintf(testFile, "%d ", stateTableAllocate[i].candidatePointNum);
	}
	//交给GPU进行并行查询
	//先传递stateTable
	//timer.stop();
	//cout << "Time 1:" << timer.elapse() << "ms" << endl;

	//timer.start();
	CUDA_CALL(cudaMemcpyAsync(this->stateTableGPU[device_idx], stateTableAllocate, sizeof(RangeQueryStateTable)*this->stateTableLength[device_idx],
		cudaMemcpyHostToDevice, stream));
	//传递完成，开始调用kernel查询
	uint8_t* resultsReturned = (uint8_t*)malloc(sizeof(uint8_t) * (this->trajNum + 1) * rangeNum);

	//timer.stop();
	//cout << "Time 2:" << timer.elapse() << "ms" << endl;

	//timer.start();
	cudaRangeQueryTestHandler((RangeQueryStateTable*)this->stateTableGPU[device_idx], this->stateTableLength[device_idx], resultsReturned, this->trajNum + 1, rangeNum, stream);
	//ofstream fp("queryResult(GTS).txt", ios_base::out);
#ifdef CHECK_CORRECT

	for (int jobID = 0; jobID <= rangeNum - 1; jobID++)
	{
		for (int traID = 0; traID <= this->trajNum; traID++)
		{
			if (resultsReturned[jobID * (this->trajNum + 1) + traID] == 1)
			{
				ResultTable[jobID][traID] = TRUE;
			}
		}
	}
#endif
	//for (vector<uint8_t>::iterator iter = resultsReturned.begin(); iter != resultsReturned.end(); iter++) {
	//	//cout << (*iter) << endl;
	//	//printf("%d\n", *iter);
	//}
	//timer.stop();
	//cout << "Time 3:" << timer.elapse() << "ms" << endl;

	//FILE *fp = fopen("resultQuery.txt", "w+");
	//for (int i = 0; i <= stateTableLength - 1; i++) {
	//	for (int j = 0; j <= stateTableAllocate[i].candidatePointNum - 1; j++) {

	//		if ((resultsReturned[i*maxPointNum + j]) == (uint8_t)(1)) {
	//			fprintf(fp,"%d\n", stateTableAllocate[i].startIdxInAllPoints + j);
	//			fprintf(fp,"%f,%f\n", allPoints[stateTableAllocate[i].startIdxInAllPoints + j].x, allPoints[stateTableAllocate[i].startIdxInAllPoints + j].y);
	//		}

	//	}
	//}
	//查询结束，善后，清空stateTable，清空gpu等
	cudaStreamDestroy(stream);
	return 0;
}

int MortonGrid::rangeQueryBatchMultiGPU(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize)
{
	MyTimer timer;
	int device_num = 2;
	vector<thread> threads_RQ;
	int rangeNumGPU[2];
	rangeNumGPU[0] = rangeNum / 2;
	rangeNumGPU[1] = rangeNum - rangeNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = rangeNumGPU[0];
	void* allocatedGPUMem[2] = { NULL };
	vector<RangeQueryStateTable> stateTableRange[2];
	stateTableRange[0].resize(rangeNum * 50000);
	stateTableRange[1].resize(rangeNum * 50000);

	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)&this->baseAddrRange[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)&this->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		allocatedGPUMem[device_idx] = this->baseAddrRange[device_idx];
		threads_RQ.push_back(thread(std::mem_fn(&MortonGrid::rangeQueryBatchGPU), this, &bounds[startIdx[device_idx
		]], rangeNumGPU[device_idx], &ResultTable[startIdx[device_idx]], resultSetSize, &stateTableRange[device_idx][0], device_idx));
	}
	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Dual GPU Time:" << timer.elapse() << "ms" << endl;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
		CUDA_CALL(cudaFree(this->stateTableGPU[device_idx]));
	}
	return 0;
}

int MortonGrid::findMatchNodeInQuadTreeGPU(MortonNode node, MBB & bound, std::vector<MortonNode>* cells, cudaStream_t stream, int queryID, int device_idx)
{
	int totalLevel = int(log2(this->cellnum) / log2(4));
	int idxInBitmap = findIdxInBitmap(node.level, node.nid);
	if (this->isLeaf[idxInBitmap])
	{
		// printf("level:%d,node:%d\n",node.level,node.nid);
		int startCellID = node.nid * int(pow(4, (totalLevel - node.level)));
		// printf("startcellID:%d\n", startCellID);
		int startIdx = this->cellPtr[startCellID].pointRangeStart;
		int pointNum = 0;
		for (int i = node.nid * int(pow(4, (totalLevel - node.level))+0.5); i <= (node.nid + 1) * int(pow(4, (totalLevel - node.level))+0.5) - 1; i++)
		{
			pointNum += this->cellPtr[i].totalPointNum;
		}
		SPoint* dataPtr = NULL;
		//如果gpu内存中没有该node的信息
		if (this->nodeAddrTable[device_idx].find(startCellID) == this->nodeAddrTable[device_idx].end())
		{
			CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum * sizeof(SPoint), cudaMemcpyHostToDevice, stream));
			dataPtr = (SPoint*)this->baseAddrRange[device_idx];
			this->nodeAddrTable[device_idx].insert(pair<int, void*>(startCellID, this->baseAddrRange[device_idx]));
			this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
		}
		//如果有，不再复制，直接用
		else
		{
			//this->stateTableRange[device_idx]->ptr = this->nodeAddrTable[device_idx].find(startCellID)->second;
			dataPtr = (SPoint*)this->nodeAddrTable[device_idx].find(startCellID)->second;
		}
		int pointsInState = 0;
		for (int idx = 0; idx < pointNum; idx += MAXPOINTINNODE) {
			this->stateTableRange[device_idx]->ptr = dataPtr;
			this->stateTableRange[device_idx]->xmin = bound.xmin;
			this->stateTableRange[device_idx]->xmax = bound.xmax;
			this->stateTableRange[device_idx]->ymin = bound.ymin;
			this->stateTableRange[device_idx]->ymax = bound.ymax;
			if (idx + MAXPOINTINNODE >= pointNum)
				pointsInState = pointNum - idx;
			else
				pointsInState = MAXPOINTINNODE;
			this->stateTableRange[device_idx]->candidatePointNum = pointsInState;
			this->stateTableRange[device_idx]->startIdxInAllPoints = startIdx + idx;
			this->stateTableRange[device_idx]->queryID = queryID;
			this->stateTableRange[device_idx] = this->stateTableRange[device_idx] + 1;
			this->stateTableLength[device_idx] = this->stateTableLength[device_idx] + 1;
			this->testCnt++;
			//printf("%d ", this->testCnt);
			dataPtr += pointsInState;
		}
	}
	else
	{
		MBB UL = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 0);
		if (bound.intersect(UL)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 0;
			findMatchNodeInQuadTreeGPU(tempNode, bound, cells, stream, queryID, device_idx);
		}
		MBB UR = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 1);
		if (bound.intersect(UR)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 1;
			findMatchNodeInQuadTreeGPU(tempNode, bound, cells, stream, queryID, device_idx);
		}
		MBB DL = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 2);
		if (bound.intersect(DL)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 2;
			findMatchNodeInQuadTreeGPU(tempNode, bound, cells, stream, queryID, device_idx);
		}
		MBB DR = this->generateMBBfromNode(node.level + 1, ((node.nid) << 2) + 3);
		if (bound.intersect(DR)) {
			MortonNode tempNode;
			tempNode.level = node.level + 1;
			tempNode.nid = (node.nid << 2) + 3;
			findMatchNodeInQuadTreeGPU(tempNode, bound, cells, stream, queryID, device_idx);
		}
	}
	return 0;
}

int MortonGrid::SimilarityQueryBatch(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue)
{
	return 0;
}

int MortonGrid::SimilarityQueryBatchCPUParallel(Trajectory * qTra, int queryTrajNum, int * EDRdistance, int kValue)
{
	return 0;
}

int MortonGrid::SimilarityMultiThreadHandler(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp>* queryQueue, Trajectory * qTra, int queryTrajNum, std::priority_queue<FDwithID, std::vector<FDwithID>, cmpBig>* EDRCalculated, int kValue, int startQueryIdx)
{
	return 0;
}

int MortonGrid::FDCalculateParallelHandeler(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp>* queue, std::map<int, int>* freqVectorQ)
{
	return 0;
}

int MortonGrid::SimilarityExecuter(SPoint * queryTra, SPoint ** candidateTra, int queryLength, int * candidateLength, int candSize, int * resultArray)
{
	return 0;
}

int MortonGrid::SimilarityQueryBatchOnGPU(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue)
{
	return 0;
}

int MortonGrid::SimilarityQueryBatchOnMultiGPU(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue)
{
	return 0;
}
