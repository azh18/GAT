#pragma once

struct CPURangeQueryResult
{
	int traid;
	float x;
	float y;
	CPURangeQueryResult *next;
};

class QueryResult
{
public:
	QueryResult();
	~QueryResult();

	CPURangeQueryResult* start;
};

