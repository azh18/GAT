#pragma once
#include <vector>

typedef std::vector<bool> CPURangeQueryResult;

typedef struct RangeQueryResultGPU {
	short jobID;
	short idx; //在candidate中第几个，以便传回去后查找
}RangeQueryResultGPU;
class QueryResult
{
public:
	QueryResult();
	~QueryResult();

	CPURangeQueryResult* start;
};

