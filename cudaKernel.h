#pragma once
#include "ConstDefine.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <assert.h>
#include <vector>

//调用CUDA，传递point地址、point数量、查询的mbb、结果地址、结果数量

#define CUDA_CALL(x) { const cudaError_t a = (x); if (a!= cudaSuccess) { printf("\nCUDA Error: %s(err_num=%d)\n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}

#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

#ifdef _CELL_BASED_STORAGE
int cudaRangeQueryHandler(int* candidateCells, int* rangeStarts, int* rangeEnds, int candidateCellNum, float xmin, float ymin, float xmax, float ymax, Point*& resultsGPU, int& resultNum, Point *pointsPtrGPU, Point *&result);
int putCellDataSetIntoGPU(Point* pointsPtr, Point*& pointsPtrGPU, int pointNum);
#endif

int handleEDRdistance(SPoint *queryTra, SPoint **candidateTra, int candidateNum, int queryLength, int *candidateLength, int *result);
void CUDAwarmUp();
int cudaRangeQueryTestHandler(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxQueryNum, cudaStream_t stream);

int EDRDistance_Batch_Handler(int queryTaskNum, TaskInfoTableForSimilarity* taskInfoTable, SPoint *queryTra, int* queryTraOffset, SPoint** candidateTraOffsets, int* queryLength, int *candidateLength, int *result, cudaStream_t *stream);

int Similarity_Pruning_Handler(intPair* queryFVGPU, intPair* FVinfo, intPair* FVTable, intPair* SubbedArray, intPair* SubbedArrayOffset, int SubbedArrayJump, int queryCellLength, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance, cudaStream_t stream);

int cudaRangeQuerySTIGHandler(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxJobNum, cudaStream_t stream);
int cudaRangeQueryTestHandlerNonMorton(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxJobNum, cudaStream_t stream, int* offset, int* offsetLen, int* offsetNum, int* offsetInOffset);
//__global__ void testSharedMemory();
