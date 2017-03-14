#pragma once
#define MAX_TRAJ_SIZE 100000
#define MAXLENGTH 1024
//MAXGAP是最大轨迹内时间间隔，如果超过这个间隔应该被视为两条轨迹
#define MAXGAP 3600

#include <stdio.h>
#include <string>

//test:以cell为基础存储
#define _CELL_BASED_STORAGE


typedef struct Point {
	float x;
	float y;
	uint32_t time;
	uint32_t tID;
}Point;