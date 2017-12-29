// Data_Structure_Test.cpp : ¶¨Òå¿ØÖÆÌ¨Ó¦ÓÃ³ÌÐòµÄÈë¿Úµã¡£
//

#ifndef WIN32
#include <unistd.h>
#else
#include <Windows.h>
#define sleep(x) Sleep(x*1000)
#endif

#include <iostream>
#include <fstream>
#include <map>
#include "CalcInAxis.h"
#include "PreProcess.h"
#include "Trajectory.h"
#include <vector>
#include "Grid.h"
#include "cudaKernel.h"
#include "SystemTest.h"
#include "STIG.h"
#include "MortonGrid.h"


using namespace std;



map<string, tidLinkTable*> vidTotid;
map<string, tidLinkTable*>::iterator iter;

//global
Trajectory* tradb;
string baseDate = "2014-07-01";
void* baseAddrGPU = NULL;

Trajectory splitTrajectory(const Trajectory* cand, int parts)
{
	Trajectory t;
	int candLen = cand->length;
	t.tid = 888888;
	t.length = 0;
	for (int i = 0; i <= candLen - 2;i++)
	{
		SamplePoint p1 = cand->points[i];
		SamplePoint p2 = cand->points[i + 1];
		switch (parts)
		{
		case 1:
			t.points.push_back(p1);
			t.length += 1;
			break;
		case 2:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.5, p1.lat + (p2.lat - p1.lat)*0.5, 1000, 888888));
			t.length += 2;
			break;
		case 3:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.33, p1.lat + (p2.lat - p1.lat)*0.33, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.66, p1.lat + (p2.lat - p1.lat)*0.66, 1000, 888888));
			t.length += 3;
			break;
		case 4:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.25, p1.lat + (p2.lat - p1.lat)*0.25, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.5, p1.lat + (p2.lat - p1.lat)*0.5, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.75, p1.lat + (p2.lat - p1.lat)*0.75, 1000, 888888));
			t.length += 4;
			break;
		case 5:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.2, p1.lat + (p2.lat - p1.lat)*0.2, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.4, p1.lat + (p2.lat - p1.lat)*0.4, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.6, p1.lat + (p2.lat - p1.lat)*0.6, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.8, p1.lat + (p2.lat - p1.lat)*0.8, 1000, 888888));
			t.length += 5;
			break;
		default:
			throw("error number of part");
			break;
		}
	}
	return t;
}

Trajectory reduceTrajectory(const Trajectory* cand, int parts)
{
	Trajectory t;
	int candLen = cand->length;
	t.tid = 888888;
	t.length = 0;
	for (int i = 0; i <= candLen - 1; i+=parts)
	{
		t.points.push_back(cand->points[i]);
		t.length++;
	}
	return t;
}

int main()
{
	int WriteTrajectoryToFile(string outFileName, int numTra);
	cout << "Hello world!" << endl;
	cout << sizeof(Cell) << endl;
	//ifstream fin;
	//float lon1, lat1, lon2, lat2;
	//lat1 = +34.15265;
	//lon1 = +113.10221;
	//lat2 = +35.15221;
	//lon2 = +113.10222;
	//cout << calculateDistance(lat1, lon1, lat2, lon2) << endl;
	//zero-copy ÉùÃ÷
	//CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));

	tradb = new Trajectory[MAX_TRAJ_SIZE];
	//PreProcess pp("data_SSmall_SH.txt", "dataout.txt");
	//PreProcess pp("SH_0.txt", "dataout.txt");
	//pp.writeTraDataToFile("SH_0_OUT.txt");
	PreProcess pp;
	pp.readTraFromFormatedFile("SH_0_OUT.txt");
	sleep(1);
	//cout << WriteTrajectoryToFile("dataOut.txt", pp.maxTid) << endl;
	cout << "read trajectory success!" << endl << "Start building cell index" << endl;
	//for (int i = 1; i <= 10000;i++)
	//{
	//	printf("%d,%d\t", i, tradb[i].length);
	//}
	printf("cellCV=2\n");
	Grid* g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 2);
	g->addDatasetToGrid(tradb, pp.maxTid);
	//Grid *g;
	//sleep(1);
	STIG *stig = new STIG();
	stig->initial(10240, 2, tradb, pp.maxTid);
	//STIG *stig;
	//sleep(1);
	FSG *fsg = new FSG(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN);
	fsg->addDatasetToGrid(tradb, pp.maxTid);
	
	MortonGrid *mgrid = new MortonGrid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 2);
	mgrid->addDatasetToGrid(tradb, pp.maxTid);
	
	//sleep(1);
	//int* EDRTable = new int[pp.maxTid*fsg->cellnum];
	//for(int i=0;i<=pp.maxTid*fsg->cellnum-1;i++)
	//printf("%d",pp.maxTid*fsg->cellnum*sizeof(int)/1024/1024);
	//FSG *fsg;
	//sleep(1);
	//return 0;
	cout << "build cell index success!" << endl;
	//delete[] tradb;
	int count = 0;
	//for (int i = 0; i <= g->cellnum - 1; i++)
	//{
	//	if (g->cellPtr[i].subTraNum == 0)
	//		count++;
	//}
	//cout << "zero num:" << count << "total" << g->cellnum << endl;
	//int temp[7] = { 553,554,555,556,557,558,559 };
	//int sizetemp = 7;
	//g->writeCellsToFile(temp, sizetemp, "111.txt");
	
	SystemTest test(tradb, g, stig, fsg, mgrid);
	// test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228 ), 80);
	// test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);
	// test.rangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228 ), 80);
	// test.rangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228 ), 80);
	// test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228 ), 80);
	// test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);
	// test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228 ), 80);
	// test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228 ), 80);
	// test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228 ), 80);
	// test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);
	// test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228 ), 80);
	// test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228 ), 80);

	//  test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 20);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
	// test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 140);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 20);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
	// test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
	// test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 140);
	// test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 20);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
	// test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 140);
	//test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	//test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	//test.similarityQueryTest(47, 20, 5);
	//
	/*
	test.similarityQueryTest(19358, 10, 5);
	test.similarityQueryTest(1299, 10, 5);
	test.similarityQueryTest(33708, 10, 5);
	test.similarityQueryTest(67630, 10, 5);
	test.similarityQueryTest(99263, 10, 5);
	test.similarityQueryTest(1309, 10, 5);
	test.similarityQueryTest(10702, 10, 5);
	test.similarityQueryTest(74853, 10, 5);
	test.similarityQueryTest(104728, 10, 5);
	test.similarityQueryTest(149443, 10, 5);
	test.similarityQueryTest(66375, 10, 5);	
	test.similarityQueryTest(149273, 10, 5);
	test.similarityQueryTest(200276, 10, 5);
	test.similarityQueryTest(229797, 10, 5);
	test.similarityQueryTest(249730, 10, 5);
	test.similarityQueryTest(21032, 10, 5);
	test.similarityQueryTest(51991, 10, 5);
	test.similarityQueryTest(69468, 10, 5);
	test.similarityQueryTest(92757, 10, 5);
	*/
	//test.similarityQueryTest(splitTrajectory(24269, 4), 20, 15);
	printf("test on different length:--------------------------\n");
	// test.similarityQueryTest(splitTrajectory(24269, 1), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 2), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 3), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 15);
	//int TRAID = 32519;
	int TRAID = 84199;
	//test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 1), 40, 15);
	Trajectory d = reduceTrajectory(&tradb[TRAID], 4);
	//Trajectory d = tradb[TRAID];
	
	//test.similarityQueryTest(splitTrajectory(&d, 1), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 2), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 3), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 4), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);

	// test.similarityQueryTest(splitTrajectory(&d, 1), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 2), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 3), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 4), 40, 15);
	//test.similarityQueryTest(splitTrajectory(&d, 3), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 2), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 4), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 6), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 8), 40, 15);
	printf("test on performance:--------------------------\n");
	test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);

	test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 80);

	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);

	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);

	test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);

	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 10, 5);
	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 20, 5);
	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 30, 5);
	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 40, 5);
	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 50, 5);
	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 60, 5);
	// test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 70, 5);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 20, 25);
	test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 60, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 80, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 100, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 20, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 30, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 40, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 50, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 60, 15);
	
	// // test.similarityQueryTest(21032, 30, 5);
	// // test.similarityQueryTest(21032, 40, 5);
	// // test.similarityQueryTest(21032, 50, 5);
	// // test.similarityQueryTest(21032, 60, 5);
	//printf("kkkkkkkkkkkkkkkkkkkkkkkkkkkkkk\n");
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 5);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 10);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 15);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 20);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 5);
	// 
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 10);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 20);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 25);
	/*
	test.similarityQueryTest(47, 10, 5);
	test.similarityQueryTest(47, 20, 5);
	test.similarityQueryTest(47, 30, 5);
	test.similarityQueryTest(47, 40, 5);
	test.similarityQueryTest(47, 50, 5);
	test.similarityQueryTest(47, 60, 5);
	test.similarityQueryTest(47, 70, 5);
	test.similarityQueryTest(47, 20, 5);
	test.similarityQueryTest(47, 20, 10);
	test.similarityQueryTest(47, 20, 15);
	test.similarityQueryTest(47, 20, 20);
	test.similarityQueryTest(47, 20, 25);
	test.similarityQueryTest(47, 20, 30);
	test.similarityQueryTest(47, 20, 35);
	*/


	//delete g;
	printf("test on different CV:--------------------------\n");
	// printf("cellCV=0\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 0);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test1(tradb, g, stig, fsg);
	// test1.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=2\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 2);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test2(tradb, g, stig, fsg);
	// test2.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=4\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 4);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test3(tradb, g, stig, fsg);
	// test3.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=6\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 6);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test4(tradb, g, stig, fsg);
	// test4.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=8\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 8);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test5(tradb, g, stig, fsg);
	// test5.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;

	printf("Finished.\n");
	//CPURangeQueryResult* resultTable = NULL;
	//int RangeQueryResultSize = 0;
	//MBB mbbArray[1000];
	//int* resultSize = NULL;
	//for (int i = 0; i <= 999; i++)
	//	mbbArray[i] = MBB(121.1, 31.1, 121.3, 31.3);
	//MyTimer timer;
	//timer.start();
	//g->rangeQueryBatch(mbbArray, 1000, resultTable, resultSize);
	//timer.stop();
	//cout << "CPU Time:" << timer.elapse() << "ms" << endl;

	//CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), 512 * 1024 * 1024));
	//void* baseAddr = baseAddrGPU;
	//timer.start();
	//g->rangeQueryBatchGPU(mbbArray, 1000, resultTable, resultSize);
	//timer.stop();
	//cout << "GPU Time:" << timer.elapse() << "ms" << endl;
	//CUDA_CALL(cudaFree(baseAddr));
	//baseAddrGPU = NULL;
	//Trajectory* qTra = new Trajectory[100];
	//for (int i = 0; i <= 99; i++)
	//{
	//	qTra[i] = tradb[47]; // length is 1024
	//}
	////for (int i = 1; i <= 9999;i++)
	////{
	////	if (tradb[i].length > 600)
	////		printf("tra:%d,length:%d\n", i, tradb[i].length);
	////}

	////Similarity on CPU
	//int* simiResult = new int[10 * 100];
	//g->SimilarityQueryBatch(qTra, 2, simiResult,50);
	//for (int i = 0; i <= 1; i++) {
	//	cout << "Trajectory:" << i << endl;
	//	for (int j = 40; j <= 49; j++) {
	//		cout << simiResult[i * 50 + j] << "\t" << endl;
	//	}
	//}
	//delete[] simiResult;



	////Similarity on GPU
	//simiResult = new int[10 * 100];
	//g->SimilarityQueryBatchOnGPU(qTra, 2, simiResult, 50);
	//for (int i = 0; i <= 1; i++)
	//{
	//	cout << "Trajectory:" << i << endl;
	//	for (int j = 40; j <= 49; j++)
	//	{
	//		cout << simiResult[i * 50 + j] << "\t" << endl;
	//	}
	//}


	//g->rangeQuery(MBB(121.4, 31.15, 121.6, 31.25), resultTable, &RangeQueryResultSize);
	//g->rangeQueryGPU(MBB(121.4, 31.15, 121.6, 31.25), resultTable, &RangeQueryResultSize);

	/*
	²âÊÔEDR DistanceµÄGPU°æ±¾
	17.3.28
	*/
	//Trajectory **testTra = (Trajectory**)malloc(sizeof(Trajectory*) * 5000);
	//for (int i = 2; i <= 5001; i++) {
	//	testTra[i-2] = &tradb[i];
	//}
	//float *EDRdistance = (float*)malloc(sizeof(float) * 5000);
	//g->SimilarityQuery(tradb[2], testTra, 5000, EDRdistance);
	/*
	²âÊÔ
	*/


	//ofstream ftest;




	getchar();
	getchar();
	getchar();

	return 0;
}

int WriteTrajectoryToFile(string outFileName, int numTra)
{
	ofstream fout;
	fout.open(outFileName, ios_base::out);
	for (int i = 1; i <= numTra; i++)
	{
		fout << i << ": ";
		for (int j = 0; j <= tradb[i].length - 1; j++)
		{
			fout << tradb[i].points[j].lon << "," << tradb[i].points[j].lat << ";";
		}
		fout << endl;
	}
	fout.close();
	return 1;
}
