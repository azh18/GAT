#include "SystemTest.h"



SystemTest::SystemTest()
{
}


SystemTest::~SystemTest()
{
}

SystemTest::SystemTest(Trajectory* tradb, Grid* g, STIG *stig)
{
	this->tradb = tradb;
	this->g = g;
	this->stig = stig;
}

int SystemTest::rangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	CPURangeQueryResult* resultTable = NULL;
	MBB mbbArray[5000];
	int* resultSize = NULL;
	for (int i = 0; i <= 4999; i++)
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;
	timer.start();
	g->rangeQueryBatch(mbbArray, rangeQueryNum, resultTable, resultSize);
	timer.stop();
	cout << "CPU Time:" << timer.elapse() << "ms" << endl;


	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), 512 * 1024 * 1024));
	void* baseAddr = baseAddrGPU;
	timer.start();
	g->rangeQueryBatchGPU(mbbArray, rangeQueryNum, resultTable, resultSize);
	timer.stop();
	cout << "GPU Time:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(baseAddr));
	return 0;
}

int SystemTest::similarityQueryTest(int similarityScale, int similarityKValue)
{
	baseAddrGPU = NULL;
	Trajectory* qTra = new Trajectory[similarityScale];
	for (int i = 0; i <= similarityScale-1; i++)
	{
		qTra[i] = tradb[47]; // length is 1024
	}
	//for (int i = 1; i <= 9999;i++)
	//{
	//	if (tradb[i].length > 600)
	//		printf("tra:%d,length:%d\n", i, tradb[i].length);
	//}

	//Similarity on CPU
	int* simiResult = new int[similarityKValue * similarityScale];
	//g->SimilarityQueryBatch(qTra, similarityScale, simiResult, similarityKValue);
	printf("multi-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
	g->SimilarityQueryBatchCPUParallel(qTra, similarityScale, simiResult, similarityKValue);
	/*
	for (int i = 0; i <= similarityScale-1; i++) {
		cout << "Trajectory:" << i << endl;
		for (int j = 0; j <= similarityKValue-1; j++) {
			cout << simiResult[i * similarityKValue + j] << "\t" << endl;
		}
	}
	*/
	delete[] simiResult;



	//Similarity on GPU
	simiResult = new int[similarityKValue * similarityScale];
	g->SimilarityQueryBatchOnGPU(qTra, similarityScale, simiResult, similarityKValue);
	/*
	for (int i = 0; i <= similarityScale-1; i++)
	{
		cout << "Trajectory:" << i << endl;
		for (int j = 0; j <= similarityKValue-1; j++)
		{
			cout << simiResult[i * similarityKValue + j] << "\t" << endl;
		}
	}
	*/
	delete[] simiResult;
	return 0;
}

int SystemTest::STIGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)
{
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), 512 * 1024 * 1024));
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	CPURangeQueryResult* resultTable = NULL;
	MBB mbbArray[5000];
	int* resultSize = NULL;
	for (int i = 0; i <= 4999; i++)
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;
	timer.start();
	stig->rangeQueryGPU(mbbArray, rangeQueryNum, resultTable, resultSize);
	timer.stop();
	cout << "CPU Time of STIG:" << timer.elapse() << "ms" << endl;
	return 0;
}
