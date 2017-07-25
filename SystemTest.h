#pragma once

#include "Trajectory.h"
#include "Grid.h"
#include "MBB.h"

extern void* baseAddrGPU;


class SystemTest
{
public:
	SystemTest();
	~SystemTest();
	Trajectory* tradb;
	Grid* g;
	MBB rangeQueryMBB;
	int rangeQueryNum;
	int similarityScale;
	int similarityKValue;
	SystemTest(Trajectory* tradb, Grid* g);
	int rangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum);
	int similarityQueryTest(int similarityScale, int similarityKValue);

};

