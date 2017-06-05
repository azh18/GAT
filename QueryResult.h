#pragma once

struct CPURangeQueryResult
{
	int traid;
	float x;
	float y;
	CPURangeQueryResult *next;
};

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

