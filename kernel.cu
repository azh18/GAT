
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cudaKernel.h"
#include "thrust\device_ptr.h"
#include "thrust\remove.h"
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include "WinTimer.h"
#include "ConstDefine.h"



#define CUDA_CALL(x) { const cudaError_t a = (x); if (a!= cudaSuccess) { printf("\nCUDA Error: %s(err_num=%d)\n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}

#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

inline void __getLastCudaError(const char *errorMessage, const char *file, const int line)
{
	cudaError_t err = cudaGetLastError();

	if (cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
			file, line, errorMessage, (int)err, cudaGetErrorString(err));
			exit(EXIT_FAILURE);
	}
}


//using namespace thrust;
static const int THREAD_N = 256; //每个block线程数

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);



#ifdef _CELL_BASED_STORAGE
int putCellDataSetIntoGPU(Point* pointsPtr, Point*& pointsPtrGPU, int pointNum) {
	CUDA_CALL(cudaSetDevice(0));
	CUDA_CALL(cudaMalloc((void**)&pointsPtrGPU, pointNum * sizeof(Point))); //分配数据的内存
	//debug
	//std::cout << pointNum << std::endl;
	//debug
	CUDA_CALL(cudaMemcpy(pointsPtrGPU, pointsPtr, pointNum * sizeof(Point), cudaMemcpyHostToDevice));//数据拷贝到gpu里
	return 0;
}
__global__ void cudaRangeQuery(int* rangeStarts, int* rangeEnds, int candidateCellNum, const Point* pointsPtr, const float xmin, const float ymin, const float xmax, const float ymax, const int *resultOffset, Point* resultPtrCuda) {
	int cellNo = blockIdx.x; //candidate里面第几个cell 0,1,2,....
	if (cellNo >= candidateCellNum) return;
	int tid = threadIdx.x;
	if (tid >= 256) return;
	int pointNum = rangeEnds[cellNo] - rangeStarts[cellNo] + 1;//block要处理的这个cell有这么多个点
	const int offset = rangeStarts[cellNo];
	for (int i = tid; i <= pointNum - 1; i += THREAD_N) {
		float x = pointsPtr[offset + i].x;
		float y = pointsPtr[offset + i].y;
		uint32_t tid = pointsPtr[offset + i].tID;
		uint32_t time = pointsPtr[offset + i].time;
		if (x <= xmax &&x >= xmin&&y <= ymax&&y >= ymin) {
			resultPtrCuda[resultOffset[cellNo] + i].x = x;
			resultPtrCuda[resultOffset[cellNo] + i].y = y;
			resultPtrCuda[resultOffset[cellNo] + i].tID = tid;
			resultPtrCuda[resultOffset[cellNo] + i].time = time;
		}
		else
			resultPtrCuda[resultOffset[cellNo] + i].tID = -1;
	}
}
int cudaRangeQueryHandler(int* candidateCells, int* rangeStarts, int* rangeEnds, int candidateCellNum,float xmin, float ymin, float xmax, float ymax, Point*& resultsGPU, int& resultNum,Point *pointsPtrGPU,Point *&result) {
	//第一个参数暂时没有用，注意这里candidatecells[i]已经不再是cell的id，由于只存非空
	//第四个参数表示非空的cell个数
	//计算candidate集中的点的总数，gpu内开辟相同大小的空间做flag。rangestart和rangeend是相应candidatecell内的采样点在AllPoints的起始下标和终止下标
	//倒数第三个和倒数第二个参数是输出，一个是保存结果的GPU地址，第二个是结果的个数
	//PointsPtrGPU是数据集在gpu的地址
	//MyTimer timer1;
	//timer1.start();
	int counter = 0;
	int *resultOffset = (int*)malloc(sizeof(int)*candidateCellNum);
	//std::cout << candidateCellNum << ":"<<std::endl;
	for (int i = 0; i <= candidateCellNum - 1; i++) {
		resultOffset[i] = counter;
		////debug
		//std::cout << "(" << rangeStarts[i] << "," << rangeEnds[i] << ");"<<"["<<resultOffset[i]<<"]";
		////debug
		counter += rangeEnds[i] - rangeStarts[i] + 1;
	}
	int totalPointNumInCandidate = counter;


	int* candidateCellsCuda = NULL, *rangeStartsCuda = NULL, *rangeEndsCuda = NULL, *resultOffsetCuda = NULL;

	CUDA_CALL(cudaMalloc((void**)&resultsGPU, sizeof(Point)*totalPointNumInCandidate));
	//将range和cell信息写入gpu
	//CUDA_CALL(cudaMalloc((void**)&candidateCellsCuda, sizeof(int)*candidateCellNum));
	CUDA_CALL(cudaMalloc((void**)&rangeStartsCuda, candidateCellNum*sizeof(int)));
	//std::cout << "\n" << candidateCellNum*sizeof(int) << "\n";
	CUDA_CALL(cudaMalloc((void**)&rangeEndsCuda, candidateCellNum*sizeof(int)));
	CUDA_CALL(cudaMalloc((void**)&resultOffsetCuda, candidateCellNum*sizeof(int)));
	//CUDA_CALL(cudaMemcpy(candidateCellsCuda, candidateCells, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(rangeStartsCuda, rangeStarts, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMemcpy(rangeEndsCuda, rangeEnds, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(resultOffsetCuda, resultOffset, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	////debug
	//CUDA_CALL(cudaMemcpy( rangeStarts, rangeStartsCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaMemcpy( rangeEnds, rangeEndsCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaMemcpy( resultOffset, resultOffsetCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//for (int i = 0; i <= candidateCellNum - 1; i++) {
	//	//debug
	//	std::cout << "(" << rangeStarts[i] << "," << rangeEnds[i] << ");" << "[" << resultOffset[i] << "]";
	//	//debug
	//}
	//debug
	//timer1.stop();
	//std::cout << timer1.ticks() << std::endl;
	//timer1.start();
	//调用kernel，如果某点满足条件，在相应位置写入在AllPoints中的下标，否则写入-1
	//每个cell分配给一个block
	cudaRangeQuery <<<candidateCellNum, THREAD_N >>>(rangeStartsCuda, rangeEndsCuda, candidateCellNum, pointsPtrGPU, xmin, ymin, xmax, ymax, resultOffsetCuda, resultsGPU);
	//kernel调用结束，结果保存在idxsGPU中，如果符合条件，对应元素内容是在AllPoints的下标，如果不符合，内容为-1
	//CUDA_CALL(cudaFree(candidateCellsCuda));
	CUDA_CALL(cudaFree(rangeStartsCuda));
	CUDA_CALL(cudaFree(rangeEndsCuda));
	CUDA_CALL(cudaFree(resultOffsetCuda));
	//getLastCudaError("Error in Calling 'kernel'");
	//使用Thrust删除所有-1，得到最终结果
	//timer1.stop();
	//std::cout << timer1.ticks() << std::endl;

	//串行结果合并
	//test

	//timer1.start();
	Point *resultset = NULL;
	resultset = (Point*)malloc(totalPointNumInCandidate*sizeof(Point));
	CUDA_CALL(cudaMemcpy(resultset, resultsGPU, sizeof(Point)*totalPointNumInCandidate, cudaMemcpyDeviceToHost));
	std::vector<Point> *resultPoint = new std::vector<Point>;
	for (int i = 0; i <= totalPointNumInCandidate - 1; i++) {
		if (resultset[i].tID != -1)
		{
			resultPoint->push_back(resultset[i]);
		}
	}
	result = &resultPoint->at(0);
	free(resultset);
	//test
	//timer1.stop();
	//std::cout << timer1.ticks() << std::endl;
	
	//并行结果合并
	
	//thrust::device_ptr<int> idxsPtr = thrust::device_pointer_cast(idxsGPU);
	//int a;
	//cudaMemcpy(&a, idxsGPU, 1, cudaMemcpyDeviceToHost);
	//size_t num = thrust::remove(idxsPtr, idxsPtr + totalPointNumInCandidate-1, -1) - idxsPtr;
	//int *result = (int*)malloc(sizeof(int)*num);
	//thrust::copy(idxsPtr, idxsPtr + num, result);
	//resultNum = num;
	//resultIdx = result;

	//CUDA_CALL(cudaFree(idxsGPU));


	return 0;
}





#else
int cudaRangeQueryHandler(Point* pointsPtr, int pointNum, float xmin, float ymin, float xmax, float ymax, Point*& resultsPtr, int& resultNum) {
	Point* pointsPtrCuda = NULL;
	Point* resultPtrCuda = NULL;
	CUDA_CALL(cudaMalloc((void**)&pointsPtrCuda, pointNum * sizeof(Point))); //分配数据的内存
	CUDA_CALL(cudaMalloc((void**)&resultPtrCuda, pointNum * sizeof(Point))); //gpu内存储结果的地方
	CUDA_CALL(cudaMemcpy(pointsPtrCuda, pointsPtr, pointNum * sizeof(Point), cudaMemcpyHostToDevice));//数据拷贝到gpu里

																									  //调用核函数处理数据，结果放在gpu中

																									  //取回数据，处理
	return 0;
}
#endif


__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

//int main()
//{
//    const int arraySize = 5;
//    const int a[arraySize] = { 1, 2, 3, 4, 5 };
//    const int b[arraySize] = { 10, 20, 30, 40, 50 };
//    int c[arraySize] = { 0 };
//
//    // Add vectors in parallel.
//    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "addWithCuda failed!");
//        return 1;
//    }
//
//    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
//        c[0], c[1], c[2], c[3], c[4]);
//
//    // cudaDeviceReset must be called before exiting in order for profiling and
//    // tracing tools such as Nsight and Visual Profiler to show complete traces.
//    cudaStatus = cudaDeviceReset();
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaDeviceReset failed!");
//        return 1;
//    }
//
//    return 0;
//}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
