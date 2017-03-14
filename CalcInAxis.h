#ifndef CALCINAXIS_H_INCLUDED
#define CALCINAXIS_H_INCLUDED

#include "math.h"
#include "ConstDefine.h"

//http://www.cnblogs.com/zhoug2020/p/3950933.html:距离算法
float calculateDistance(float LatA,float LonA,float LatB,float LonB)
{
    float Pi = 3.1415926;
    float R = 6371004;//单位米
    float MLatA,MLatB,MLonA,MLonB;
    MLatA = 90 - LatA;
    MLatB = 90 - LatB;
    MLonA = LonA;
    MLonB = LonB;
    float C = sin(MLatA)*sin(MLatB)*cos(MLonA - MLonB) + cos(MLatA)*cos(MLatB);
    float Distance = R*acos(C)*Pi/180;
    return Distance;
}

#endif // CALCINAXIS_H_INCLUDED
