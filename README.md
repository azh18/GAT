# GAT
A GPU-accelerated framework for processing trajectory queries, supporting two basic kinds of queries: range query and top-k similarity query.

## prerequisits ##
CUDA 8.0
GCC 4.8.5 with c++11 support

## Compile ##
VS2015: open the project file

Linux: `make all -jN`

where N is the number of cores of your CPU

## Configuration ##
open multi-GPU mode: `#define USE_MULTIGPU` in ConstDefine.h

parameters: in ConstDefine.h

## Query Tasks Setting ##

invoke following functions in main.cpp:

**range query in GT-quadtree:** `rangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)`

**range query in STIG:** `STIGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)`

**range query in FSG:** `FSGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)`

**top-k similarity in GT-quadtree on both CPU&GPU and EDR on CPU:** `similarityQueryTest(int queryTrajNo, int similarityScale, int similarityKValue)`


## Run ##
Just directly execute the compiled and linked executive file.

## Process ##
17/8/2: testing whether STIG works on multi-GPU environment

17/8/7: add multi-thread CPU version of range query

17/9/4: get wanted results, use const range query configuration to make sure the equality

17/10/9: finish the paper!
