#include "SamplePoint.h"

SamplePoint::SamplePoint()
{
    //ctor
    lon = 0;
    lat = 0;
    time = 0;
    tid = 0;
}

SamplePoint::SamplePoint(float lon1,float lat1,int time1,int tid1)
{
    //ctor
    this->lat = lat1;
    this->lon = lon1;
    this->time = time1;
    this->tid = tid1;
}

SamplePoint::~SamplePoint()
{
    //dtor
}
