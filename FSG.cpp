#include "FSG.h"
#include "cudaKernel.h"
#include "Grid.h"

extern Trajectory* tradb;

FSG::FSG(const MBB& mbb, float val_cell_size)
{
	range = mbb;
	cell_size = val_cell_size;
	cell_num_x = (int)((mbb.xmax - mbb.xmin) / val_cell_size) + 1;
	cell_num_y = (int)((mbb.ymax - mbb.ymin) / val_cell_size) + 1;
	cellnum = cell_num_x*cell_num_y;
	//找到最佳的长宽

	cellPtr = new Cell[cellnum];

	//注意cell编号是从(xmin,ymax)开始的，而不是(xmin,ymin)
	for (int i = 0; i <= this->cell_num_y - 1; i++)
	{
		for (int j = 0; j <= this->cell_num_x - 1; j++)
		{
			int cell_idx = i*this->cell_num_x + j;
			this->cellPtr[cell_idx].initial(i, j, MBB(this->range.xmin + this->cell_size*j, this->range.ymax - this->cell_size*(i + 1), this->range.xmin + this->cell_size*(j + 1), this->range.ymax - this->cell_size*(i)));
			//g->cellPtr[cell_idx].initial(i, j, MBB(g->range.xmin + g->cell_size*j, g->range.ymax - g->cell_size*(i+1), g->range.xmin + g->cell_size*(j + 1), g->range.ymax - g->cell_size*(i)));
		}
	}
}

int FSG::addTrajectoryIntoCell(Trajectory & t)
{
	if (t.length == 0)
		return 1;//空轨迹
	SamplePoint p = t.points[0];
	int lastCellNo = WhichCellPointIn(p);
	int lastCellStartIdx = 0;
	int nowCellNo;
	//cell based traj生成，记得转换后free！
	std::vector<int>* tempCellBasedTraj = new std::vector<int>;
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
			}
			//否则，上一个和这个cell都要添加
			else
			{
				tempCellNum += 2;
				tempCellBasedTraj->push_back(lastCellNo);
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				cellPtr[nowCellNo].addSubTra(t.tid, i, i, 1);
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

				lastCellNo = nowCellNo;
				lastCellStartIdx = i;
			}
		}
	}
	delete tempCellBasedTraj;
	return 0;
}

int FSG::WhichCellPointIn(SamplePoint p)
{
	int row = (int)((this->range.ymax - p.lat) / this->cell_size);
	int col = (int)((p.lon - this->range.xmin) / this->cell_size);
	return row*this->cell_num_x + col;
}

int FSG::addDatasetToGrid(Trajectory * db, int traNum)
{
	this->trajNum = traNum;
	//生成frequency vector
	//注意，轨迹编号从1开始
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
			std::cerr << "Grid.cpp: something wrong in total point statistic" << std::endl;
	}
	return 0;
}

int FSG::rangeQueryBatchGPU(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize, RangeQueryStateTable * stateTableAllocate, int device_idx)
{
#ifdef CHECK_CORRECT
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		ResultTable[i].resize(this->trajNum + 1);
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
	// for each query, generate the nodes:
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	for (int j = 0; j <= rangeNum - 1; j++)
	{
		//find candidate cells
		int g1, g2, g3, g4; //box的顶点所在行列
		int a, b, c, d;//box的顶点所在格子号
		int *candidatesCellID = NULL, *resultsCellID = NULL, *directResultsCellID = NULL;//候选格子，默认为空
		int m, n;//m为grid列数，n为行数
		int candidateSize = 0;//candidate个数
		int DirectresultSize = 0;//结果个数
		int counter = 0;//计数器
		m = this->cell_num_x;
		n = this->cell_num_y;
		g1 = (int)((bounds[j].xmin - this->range.xmin) / this->cell_size);
		g2 = (int)((bounds[j].xmax - this->range.xmin) / this->cell_size);
		g3 = (int)((this->range.ymax - bounds[j].ymax) / this->cell_size);
		g4 = (int)((this->range.ymax - bounds[j].ymin) / this->cell_size);
		a = g1 + g3*m;
		b = g2 + g3*m;
		c = g1 + g4*m;
		d = g2 + g4*m;

		if (a == b)
		{
			candidateSize = (c - a) / m + 1;
		}
		else
		{
			if (a == c)
				candidateSize = (b - a) + 1;
			else
				candidateSize = ((c - a) / m + 1) * 2 + (b - a + 1) * 2 - 4;
		}

		//把bounding box的线经过的cell加入candidates
		//这些candidate是等待检查的
		candidatesCellID = (int*)malloc(sizeof(int)*candidateSize);
		counter = 0;
		for (int i = a; i <= b; i++)
		{
			candidatesCellID[counter] = i;
			counter++;
		}
		if (g4 - g3 >= 1)
		{
			for (int i = c; i <= d; i++)
			{
				candidatesCellID[counter] = i;
				counter++;
			}
		}
		if (g4 - g3 >= 2)
		{
			for (int i = a + m; i <= a + (g4 - g3 - 1)*m; i = i + m)
			{
				candidatesCellID[counter] = i;
				counter++;
			}
			if (g2 > g1)
			{
				for (int i = b + m; i <= b + (g4 - g3 - 1)*m; i = i + m)
				{
					candidatesCellID[counter] = i;
					counter++;
				}
			}
		}
		if (counter != candidateSize)
			std::cerr << "size error in range query candidates cell" << counter << "," << candidateSize << std::endl;

		//一些直接是result
		if ((b>a) && (g4>g3))
			DirectresultSize = (b - a - 1)*(g4 - g3 - 1);
		else
			DirectresultSize = 0;


		counter = 0;
		directResultsCellID = (int*)malloc(DirectresultSize * sizeof(int));
		if (b >= a + 2 && c >= a + 2 * m)
		{
			for (int i = a + 1; i <= b - 1; i++)
			{
				for (int j = 1; j <= g4 - g3 - 1; j++)
				{
					directResultsCellID[counter] = i + j*m;
					counter++;
				}
			}
		}
		if (counter != DirectresultSize)
			std::cerr << "size error in range query directresult cell" << counter << "," << DirectresultSize << std::endl;

		//复制数据，构造pair
		for (int i = 0; i <= candidateSize - 1;i++)
		{
			int startIdx = this->cellPtr[candidatesCellID[i]].pointRangeStart;
			int pointNum = this->cellPtr[candidatesCellID[i]].pointRangeEnd - this->cellPtr[candidatesCellID[i]].pointRangeStart + 1;
			if (this->nodeAddrTable[device_idx].find(candidatesCellID[i]) == this->nodeAddrTable[device_idx].end()) {
				CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum*sizeof(SPoint), cudaMemcpyHostToDevice, stream));
				this->nodeAddrTable[device_idx][candidatesCellID[i]] = this->baseAddrRange[device_idx];
				this->stateTableRange[device_idx]->ptr = this->baseAddrRange[device_idx];
				this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
			}
			else
			{
				this->stateTableRange[device_idx]->ptr = this->nodeAddrTable[device_idx][candidatesCellID[i]];
			}
			this->stateTableRange[device_idx]->xmin = bounds[j].xmin;
			this->stateTableRange[device_idx]->xmax = bounds[j].xmax;
			this->stateTableRange[device_idx]->ymin = bounds[j].ymin;
			this->stateTableRange[device_idx]->ymax = bounds[j].ymax;
			this->stateTableRange[device_idx]->candidatePointNum = pointNum;
			this->stateTableRange[device_idx]->startIdxInAllPoints = startIdx;
			this->stateTableRange[device_idx]->queryID = j;
			this->stateTableRange[device_idx] = this->stateTableRange[device_idx] + 1;
			this->stateTableLength[device_idx] = this->stateTableLength[device_idx] + 1;
		}
		
		for (int i = 0; i <= DirectresultSize - 1; i++)
		{
			int startIdx = this->cellPtr[directResultsCellID[i]].pointRangeStart;
			int pointNum = this->cellPtr[directResultsCellID[i]].pointRangeEnd - this->cellPtr[directResultsCellID[i]].pointRangeStart + 1;
			if (this->nodeAddrTable[device_idx].find(directResultsCellID[i]) == this->nodeAddrTable[device_idx].end()) {
				CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum*sizeof(SPoint), cudaMemcpyHostToDevice, stream));
				this->nodeAddrTable[device_idx][directResultsCellID[i]] = this->baseAddrRange[device_idx];
				this->stateTableRange[device_idx]->ptr = this->baseAddrRange[device_idx];
				this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
			}
			else
			{
				this->stateTableRange[device_idx]->ptr = this->nodeAddrTable[device_idx][directResultsCellID[i]];
			}
			this->stateTableRange[device_idx]->xmin = bounds[j].xmin;
			this->stateTableRange[device_idx]->xmax = bounds[j].xmax;
			this->stateTableRange[device_idx]->ymin = bounds[j].ymin;
			this->stateTableRange[device_idx]->ymax = bounds[j].ymax;
			this->stateTableRange[device_idx]->candidatePointNum = pointNum;
			this->stateTableRange[device_idx]->startIdxInAllPoints = startIdx;
			this->stateTableRange[device_idx]->queryID = j;
			this->stateTableRange[device_idx] = this->stateTableRange[device_idx] + 1;
			this->stateTableLength[device_idx] = this->stateTableLength[device_idx] + 1;
		}
		/*
		//直接是结果的，可以全部写进结果中
		for (int i = 0; i <= DirectresultSize - 1; i++)
		{
			int startIdx = this->cellPtr[directResultsCellID[i]].pointRangeStart;
			int pointNum = this->cellPtr[directResultsCellID[i]].pointRangeEnd - this->cellPtr[candidatesCellID[i]].pointRangeStart + 1;
			for (int idx = startIdx; idx < startIdx + pointNum; idx++)
				ResultTable[j][this->allPoints[idx].tID] = TRUE;
		}
		*/
	}
	//printf("StateTableLength:%d",this->stateTableLength);
	//stateTable中点的数目的最大值
	int maxPointNum = 0;
	for (int i = 0; i <= this->stateTableLength[device_idx] - 1; i++)
	{
		if (stateTableAllocate[i].candidatePointNum > maxPointNum)
			maxPointNum = stateTableAllocate[i].candidatePointNum;
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
	//for (int jobID = 0; jobID <= rangeNum - 1; jobID++)
	//{
	//	for (int traID = 0; traID <= this->trajNum; traID++)
	//	{
	//		if (resultsReturned[jobID * (this->trajNum + 1) + traID] == 1)
	//		{
	//			fp << "job " << jobID << "find" << traID << endl;
	//		}
	//	}
	//}
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

int FSG::rangeQueryBatchMultiGPU(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize)
{
	MyTimer timer;
	int device_num = 2;
	std::vector<std::thread> threads_RQ;
	int rangeNumGPU[2];
	rangeNumGPU[0] = rangeNum / 2;
	rangeNumGPU[1] = rangeNum - rangeNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = rangeNumGPU[0];
	void* allocatedGPUMem[2] = { NULL };
	std::vector<RangeQueryStateTable> stateTableRange[2];
	stateTableRange[0].resize(rangeNum * 1000);
	stateTableRange[1].resize(rangeNum * 1000);

	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)&this->baseAddrRange[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)&this->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		allocatedGPUMem[device_idx] = this->baseAddrRange[device_idx];
		threads_RQ.push_back(std::thread(std::mem_fn(&FSG::rangeQueryBatchGPU), this, &bounds[startIdx[device_idx
		]], rangeNumGPU[device_idx], &ResultTable[startIdx[1]], resultSetSize, &stateTableRange[device_idx][0], device_idx));
	}
	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	std::cout << "Dual GPU Time:" << timer.elapse() << "ms" << std::endl;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
		CUDA_CALL(cudaFree(this->stateTableGPU[device_idx]));
	}
	return 0;
}




FSG::FSG()
{
}


FSG::~FSG()
{
}
