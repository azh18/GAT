#pragma once
#include "ConstDefine.h"
#include "Trajectory.h"
#include "Grid.h"
#include "MBB.h"
#include "STIG.h"

extern void* baseAddrGPU;


class SystemTest
{
public:
	SystemTest();
	~SystemTest();
	Trajectory* tradb;
	Grid* g;
	STIG* stig;
	MBB rangeQueryMBB;
	int rangeQueryNum;
	int similarityScale;
	int similarityKValue;
	SystemTest(Trajectory* tradb, Grid* g, STIG *stig);
	int rangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum);
	int similarityQueryTest(int similarityScale, int similarityKValue);
	int STIGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum);

};

