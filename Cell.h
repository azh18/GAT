#pragma once
#include "ConstDefine.h"
#include "MBB.h"
#include <string>
#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>

using namespace std;
typedef struct subTra {
	int traID;
	short startpID;
	short endpID;
	short numOfPoint;//每个子轨迹点的个数
	subTra* next;
}subTra;

class Cell
{
public:
	Cell();
	Cell(int x, int y,const MBB& val_mbb);
	bool initial(int x, int y, const MBB& val_mbb);
	int addSubTra(int traID, int startIdx, int endIdx, int numOfPoints);
	int buildSubTraTable();//读取完所有轨迹之后用数组存储
	int writeCellToFile(string fileName);
	~Cell();


	// 这三个属性都可以计算出来
	int cell_x;//cell横坐标
	int cell_y;//cell纵坐标
	MBB mbb;

	int subTraNum; //子轨迹的个数
	int totalPointNum; //cell内点个数

	subTra subTraEntry;//建立cell过程中链表入口（这一点改掉，太浪费内存）
	subTra* subTraPtr;//当前最新数据指针位置
	subTra* subTraTable;//转化为数组后数组的入口
	ofstream fout;//文件接口

#ifdef _CELL_BASED_STORAGE
	//每个cell在point数组中数据的起始和终止位置
	int pointRangeStart;
	int pointRangeEnd;
#endif // _CELL_BASED_STORAGE

};

