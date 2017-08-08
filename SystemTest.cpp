#include "SystemTest.h"


using namespace std;

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
		rangeQueryMBB.randomGenerateMBB(mbbArray[i]);
	MyTimer timer;

	
	printf("single-core CPU range query #query=%d:\n", rangeQueryNum);
	vector<CPURangeQueryResult> rangeQueryResult;
	rangeQueryResult.resize(rangeQueryNum);
	timer.start();
	g->rangeQueryBatch(mbbArray, rangeQueryNum, &rangeQueryResult[0], resultSize);
	timer.stop();
	cout << "single-core CPU Time:" << timer.elapse() << "ms" << endl;
	

	printf("multi-core CPU range query #query=%d:\n", rangeQueryNum);
	vector<CPURangeQueryResult> rangeQueryResultMultiCPU;
	rangeQueryResultMultiCPU.resize(rangeQueryNum);
	timer.start();
	g->rangeQueryBatchMultiThread(mbbArray, rangeQueryNum, &rangeQueryResultMultiCPU[0], resultSize);
	timer.stop();
	cout << "multi-core CPU Time:" << timer.elapse() << "ms" << endl;
	

	/*
	for (int i = 0; i <= rangeQueryNum - 1;i++)
	{
		for (int traID = 1; traID <= this->g->trajNum; traID++) {
			if (rangeQueryResultMultiCPU[i][traID])
				printf("Query %d result: %d\n", i, traID);
		}
	}
	*/

	
	printf("single GPU range query #query=%d:\n", rangeQueryNum);
	CUDA_CALL(cudaSetDevice(0));

#ifdef WIN32
	CUDA_CALL(cudaMalloc((void**)(&g->baseAddrRange[0]), (long long int)512 * 1024 * 1024));
#else
	CUDA_CALL(cudaMalloc((void**)(&g->baseAddrRange[0]), (long long int)2048 * 1024 * 1024));
#endif

	void *allocatedGPUMem = g->baseAddrRange[0];
	CUDA_CALL(cudaMalloc((void**)&g->stateTableGPU[0], 512 * 1024 * 1024));
	vector<RangeQueryStateTable> stateTableRange;
	stateTableRange.resize(rangeQueryNum * 1000);
	timer.start();
	g->rangeQueryBatchGPU(mbbArray, rangeQueryNum, resultTable, resultSize, &stateTableRange[0], 0);
	timer.stop();
	cout << "Single GPU Time:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(allocatedGPUMem));
	CUDA_CALL(cudaFree(g->stateTableGPU[0]));

#ifdef USE_MULTIGPU
	printf("multi-GPU range query #query=%d:\n", rangeQueryNum);
	g->rangeQueryBatchMultiGPU(mbbArray, rangeQueryNum, resultTable, resultSize);
#else

#endif
	

	return 0;
}

int SystemTest::similarityQueryTest(int queryTrajNo, int similarityScale, int similarityKValue)
{
	baseAddrGPU = NULL;
	Trajectory* qTra = new Trajectory[similarityScale];
	// form query trajectories
	//for (int subID = 0; subID <= this->g->cellBasedTrajectory[queryTrajNo].length - 1; subID++)
	//{
	//	int idxInAllPoints = this->g->cellBasedTrajectory[queryTrajNo].startIdx[subID];
	//	memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
	//	//for (int cnt = 0; cnt <= this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID] - 1; cnt++) {
	//	//	candidateTra[i][cnt] = this->allPoints[idxInAllPoints+cnt];
	//	//}
	//	//printf("%d ", this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
	//	tempPtr += this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID];
	//}

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
	g->SimilarityQueryBatch(qTra, similarityScale, simiResult, similarityKValue);
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
	printf("one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnGPU(qTra, similarityScale, simiResult, similarityKValue);

#ifdef USE_MULTIGPU
	printf("multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPU(qTra, similarityScale, simiResult, similarityKValue);
#else

#endif
	
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
	CUDA_CALL(cudaSetDevice(0));
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	CPURangeQueryResult* resultTable = NULL;
	MBB mbbArray[5000];
	int* resultSize = NULL;
	for (int i = 0; i <= 4999; i++)
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;

	// single GPU
	printf("single GPU STIG range query #query=%d:\n", rangeQueryNum);

	void* allocatedGPUMemS = 0;
	CUDA_CALL(cudaMalloc((void**)&this->stig->baseAddrGPU[0], (long long int)512 * 1024 * 1024));
	CUDA_CALL(cudaMalloc((void**)&this->stig->stateTableGPU[0], 512 * 1024 * 1024));
	allocatedGPUMemS = this->stig->baseAddrGPU[0];
	timer.start();
	stig->rangeQueryGPU(mbbArray, rangeQueryNum, resultTable, resultSize, 0);
	timer.stop();
	cout << "single GPU Time of STIG:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(allocatedGPUMemS));
	CUDA_CALL(cudaFree(this->stig->stateTableGPU[0]));

#ifdef USE_MULTIGPU
	// multi-GPU
	int device_num = 2;
	vector<thread> threads_RQ;
	int rangeNumGPU[2];
	rangeNumGPU[0] = rangeQueryNum / 2;
	rangeNumGPU[1] = rangeQueryNum - rangeNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = rangeNumGPU[0];
	void* allocatedGPUMem[2] = { NULL };
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)&this->stig->baseAddrGPU[device_idx], (long long int)2048 * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)&this->stig->stateTableGPU[device_idx], 512 * 1024 * 1024));
		allocatedGPUMem[device_idx] = this->stig->baseAddrGPU[device_idx];
		threads_RQ.push_back(thread(std::mem_fn(&STIG::rangeQueryGPU), this->stig, &mbbArray[startIdx[device_idx
		]], rangeNumGPU[device_idx], resultTable, resultSize, device_idx));
	}
	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Dual GPU Time of STIG:" << timer.elapse() << "ms" << endl;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
		CUDA_CALL(cudaFree(this->stig->stateTableGPU[device_idx]));
	}
#else

#endif

	return 0;
}
