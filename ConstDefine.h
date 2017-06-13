#pragma once
#define MAX_TRAJ_SIZE 50000
#define MAXLENGTH 4096
//MAXGAP是最大轨迹内时间间隔，如果超过这个间隔应该被视为两条轨迹
#define MAXGAP 3600

#define EPSILON 0.001
#define MAXTHREAD 256

//每个node内包含的点的个数上限
#define MAXPOINTINNODE 1000

#include <stdio.h>
#include <string>
#include <math.h>
#include "QueryResult.h"

//test:以cell为基础存储
#define _CELL_BASED_STORAGE
//test:Similarity query based on naive grid，以定大小的grid来索引
//#define _SIMILARITY

//4+4+4+4=16bytes
typedef struct Point {
	float x;
	float y;
	uint32_t time;
	uint32_t tID;
}Point;

//4+4+4=12bytes
typedef struct SPoint {
	float x;
	float y;
	uint32_t tID;
}SPoint;

//2+2+4=8bytes
typedef struct DPoint {
	short x;
	short y;
	uint32_t tID;
}DPoint;

typedef struct cellBasedTraj {
	int *cellNo = NULL;
	int *startIdx = NULL; //相应轨迹在cell中的开始的轨迹在AllPoints中的idx
	short *numOfPointInCell = NULL;//在每个cell中对应该轨迹点的个数
	short length;
	int trajLength;
}cellBasedTraj;

typedef struct RangeQueryStateTable {
	void* ptr;
	int candidatePointNum;
	float xmin;
	float ymin;
	float xmax;
	float ymax;
	int queryID;
	int startIdxInAllPoints;
}RangeQueryStateTable;

int getIdxFromXY(int x, int y);

