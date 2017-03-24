#pragma once
#include "ConstDefine.h"
//调用CUDA，传递point地址、point数量、查询的mbb、结果地址、结果数量
#ifdef _CELL_BASED_STORAGE
int cudaRangeQueryHandler(int* candidateCells, int* rangeStarts, int* rangeEnds, int candidateCellNum, float xmin, float ymin, float xmax, float ymax, Point*& resultsGPU, int& resultNum, Point *pointsPtrGPU, Point *&result);
int putCellDataSetIntoGPU(Point* pointsPtr, Point*& pointsPtrGPU, int pointNum);
#endif