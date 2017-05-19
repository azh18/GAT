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

extern Trajectory* tradb;

class Grid
{
public:
	Grid();
	~Grid();
	Grid(const MBB& mbb, float val_cell_size);
	int addTrajectoryIntoCell(Trajectory &t);
	int WhichCellPointIn(SamplePoint p);
	int addDatasetToGrid(Trajectory* db,int traNum);
	int writeCellsToFile(int* cellNo, int cellNum,string file);
	//rangeQuery函数，输入Bounding box，输出轨迹编号和对应顺序下的采样点
	int rangeQuery(MBB & bound, CPURangeQueryResult * ResultTable, int* resultSetSize);
	int rangeQueryGPU(MBB & bound, CPURangeQueryResult * ResultTable, int* resultSetSize);
	int SimilarityQuery(Trajectory &qTra, Trajectory **candTra, int candSize, float *EDRdistance);
	static int getIdxFromXY(int x, int y);


	//Grid索引包含的坐标范围
	MBB range;
	float cell_size; //length of a cell
	int cellNum_axis; // 每行/列有几个cell
	int cellnum; //upper(area(grid)/area(cell))，保证能放下所有cell
	Cell* cellPtr; //存储cell的入口
	ofstream fout;//文件输出接口
	int totalPointNum; //grid内点个数

	vector<cellBasedTraj> cellBasedTrajectory; //cellbasedtrajectory，二元组：（cell编号数组地址，数组长度）

	Point* allPoints;//存储所有点的数组
	Point* allPointsPtrGPU;
	



};

