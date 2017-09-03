#pragma once
#include "ConstDefine.h"
#include "Trajectory.h"
#include "Grid.h"
#include "MBB.h"
#include "STIG.h"
#include "FSG.h"

extern void* baseAddrGPU;


class SystemTest
{
public:
	SystemTest();
	~SystemTest();
	Trajectory* tradb;
	Grid* g;
	STIG* stig;
	FSG* fsg;
	MBB rangeQueryMBB;
	int rangeQueryNum;
	int similarityScale;
	int similarityKValue;
	SystemTest(Trajectory* tradb, Grid* g, STIG *stig, FSG* fsg);
	int rangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum);
	int similarityQueryTest(int queryTrajNo, int similarityScale, int similarityKValue);
	int STIGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum);
	int FSGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum);

};

