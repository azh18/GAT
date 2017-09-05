#pragma once
#include "ConstDefine.h"
#include "Cell.h"
#include <map>

#include "Trajectory.h"
#include "SamplePoint.h"

class FSG
{
public:
	//Grid索引包含的坐标范围
	MBB range;
	float cell_size; //length of a cell
	int cell_num_x, cell_num_y; //横竖各有多少个cell
	int cellnum; //upper(area(grid)/area(cell))，保证能放下所有cell
	Cell* cellPtr; //存储cell的入口
	std::ofstream fout;//文件输出接口
	int totalPointNum; //grid内点个数
	int trajNum;


	SPoint* allPoints;//存储所有点的数组
	Point* allPointsPtrGPU;
	DPoint *allPointsDeltaEncoding;//Delta Encoding后的点

	//Range Query on GPU 用
	void *baseAddrRange[2];
	void *stateTableGPU[2];
	RangeQueryStateTable* stateTableRange[2];
	std::map<int, void*> nodeAddrTable[2];
	int stateTableLength[2];
	//int nodeAddrTableLength[2];

	FSG(const MBB& mbb, float val_cell_size);
	int addTrajectoryIntoCell(Trajectory &t);
	int WhichCellPointIn(SamplePoint p);
	int addDatasetToGrid(Trajectory* db, int traNum);

	//Query方法
	int rangeQueryBatchGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize, RangeQueryStateTable* stateTableAllocate, int device_idx);
	int rangeQueryBatchMultiGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);

	FSG();
	~FSG();
};

