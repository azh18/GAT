// Data_Structure_Test.cpp : 定义控制台应用程序的入口点。
//
#ifndef WIN32
#include <unistd.h>
#else
#include <Windows.h>
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


using namespace std;

map<string, tidLinkTable*> vidTotid;
map<string, tidLinkTable*>::iterator iter;

//global
Trajectory* tradb;
string baseDate = "2014-07-01";
void* baseAddrGPU = NULL;


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
	//zero-copy 声明
	//CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));

	tradb = new Trajectory[MAX_TRAJ_SIZE];
	//PreProcess pp("data_SSmall_SH.txt", "dataout.txt");
	//PreProcess pp("SH_0.txt", "dataout.txt");
	//pp.writeTraDataToFile("SH_0_OUT.txt");
	PreProcess pp;
	pp.readTraFromFormatedFile("SH_0_OUT.txt");

	//cout << WriteTrajectoryToFile("dataOut.txt", pp.maxTid) << endl;
	cout << "read trajectory success!" << endl << "Start building cell index" << endl;
	//for (int i = 1; i <= 10000;i++)
	//{
	//	printf("%d,%d\t", i, tradb[i].length);
	//}

	Grid* g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN);
	g->addDatasetToGrid(tradb, pp.maxTid);
	//Grid *g;
	//sleep(1);
	STIG *stig = new STIG();
	stig->initial(10240, 2, tradb, pp.maxTid);
	//STIG *stig;
	//sleep(1);
	FSG *fsg = new FSG(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN);
	fsg->addDatasetToGrid(tradb, pp.maxTid);
	//FSG *fsg;
	//sleep(1);
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

	SystemTest test(tradb, g, stig, fsg);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 20);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 40);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 60);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 100);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 120);
	test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 140);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 20);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 40);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 60);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 100);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 120);
	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 140);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 20);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 40);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 60);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 100);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 120);
	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 140);
	//test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	//test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	//test.similarityQueryTest(47, 20, 5);
	
	/*
	test.similarityQueryTest(47, 10, 5);
	test.similarityQueryTest(47, 20, 5);
	test.similarityQueryTest(47, 30, 5);
	test.similarityQueryTest(47, 40, 5);
	test.similarityQueryTest(47, 50, 5);
	test.similarityQueryTest(47, 60, 5);
	test.similarityQueryTest(47, 70, 5);
	*/
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
	测试EDR Distance的GPU版本
	17.3.28
	*/
	//Trajectory **testTra = (Trajectory**)malloc(sizeof(Trajectory*) * 5000);
	//for (int i = 2; i <= 5001; i++) {
	//	testTra[i-2] = &tradb[i];
	//}
	//float *EDRdistance = (float*)malloc(sizeof(float) * 5000);
	//g->SimilarityQuery(tradb[2], testTra, 5000, EDRdistance);
	/*
	测试
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
