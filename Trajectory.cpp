#include "Trajectory.h"


extern float calculateDistance(float LatA,float LonA,float LatB,float LonB);

Trajectory::Trajectory()
{

    //ctor
}

Trajectory::Trajectory(int tid,std::string vid)
{
    this->tid = tid;
    this->vid = vid;
    this->length = 0;
    //ctor
}

//往轨迹里添加采样点
//0:成功 1:超过最大数目 2:时间跨度太大，应计入下条轨迹 3:时间跨度不大，空间跨度太大，计算出来的速度大于180km/h~~50m/s，舍弃该点
int Trajectory::addSamplePoints(float lon,float lat,int time)
{
    //length指向的就是当前的idx，length超过最大值，则失败，返回1
    if(this->length>=MAXLENGTH)
    {
        return 1;
    }
    if(this->length>0)
    {
        if(time - this->points[this->length-1].time > MAXGAP)
        {
            return 2;
        }
        if((calculateDistance(lat,lon,this->points[this->length-1].lat,this->points[this->length-1].lon))/(time - this->points[this->length-1].time)>=50)
        {
            return 3;
        }
    }
    //经过检查可以加入这点到轨迹中
    this->points[this->length].lat = lat;
    this->points[this->length].lon = lon;
    this->points[this->length].tid = this->tid;
    this->points[this->length].time = time;
    this->length++;
    return 0;
}



Trajectory::~Trajectory()
{
    //dtor
}
