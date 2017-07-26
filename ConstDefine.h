#pragma once
#define MAX_TRAJ_SIZE 50000
//受到GPU可用的shared memory限制
#define MAXLENGTH 1024
//MAXGAP是最大轨迹内时间间隔，如果超过这个间隔应该被视为两条轨迹
#define MAXGAP 3600

#define EPSILON 0.001
#define MAXTHREAD 256

//每个node内包含的点的个数上限
#define MAXPOINTINNODE 1000

//在FVTable中，由于GPU显存限制，每次可以检查的FV的上限
#define N_BATCH_QUERY 2048
#define TRUE 1
#define FALSE 0
#define KSIMILARITY 80

// #define NOT_COLUMN_ORIENTED

#include <stdio.h>
#include <string>
#include <math.h>
#include "QueryResult.h"
#include <cstring>
#include <thread>


#ifdef WIN32
	#include "WinTimer.h"
#else
	#include <sys/time.h>
#endif

//test:以cell为基础存储
#define _CELL_BASED_STORAGE
//test:Similarity query based on naive grid，以定大小的grid来索引
//#define _SIMILARITY

//4+4+4+4=16bytes
typedef struct Point {
	float x;
	float y;
	int time;
	int tID;
}Point;

//4+4+4=12bytes
typedef struct SPoint {
	float x;
	float y;
	int tID;
}SPoint;

//2+2+4=8bytes
typedef struct DPoint {
	short x;
	short y;
	int tID;
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

typedef struct OffsetTable {
	int objectId;
	void *addr;
}OffsetTable;

typedef struct TaskInfoTableForSimilarity {
	int qID;
	int candTrajID;
}TaskInfoTableForSimilarity;

typedef struct intPair{
	int int_1;
	int int_2;
}intPair;

#ifdef WIN32
#else


class MyTimer
{
public:
	MyTimer(){
	};
	double iStart;
	double iEnd;

	double cpuSecond() {
		struct timeval tp;
		gettimeofday(&tp, NULL);
		return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
	}

	inline void start()
	{
		iStart = cpuSecond();
	}
	inline void stop()
	{
		iEnd = cpuSecond();
	}
	inline float elapse()
	{
		return iEnd - iStart;
	}
};
#endif

typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;

int getIdxFromXY(int x, int y);

