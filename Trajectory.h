#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include "ConstDefine.h"
#include "SamplePoint.h"


// 轨迹类，记录该轨迹的采样点等信息，最长1024
class Trajectory
{
    public:
        Trajectory();
        Trajectory(int tid,std::string vid);
        int addSamplePoints(float lon,float lat,int time);
        virtual ~Trajectory();
        int tid;
        SamplePoint points[MAXLENGTH]; //下标从0开始
        int length = 0;
        std::string vid;
        int errCounter = 0;
        SamplePoint errPointBuff[10]; //备用

    protected:

    private:
};

#endif // TRAJECTORY_H
