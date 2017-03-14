#ifndef SAMPLEPOINT_H
#define SAMPLEPOINT_H
#include "ConstDefine.h"

class SamplePoint
{
    public:
        SamplePoint();
        SamplePoint(float lon1,float lat1,int time1,int tid1);
        virtual ~SamplePoint();
        float lon;
        float lat;
        int time;
        int tid;

    protected:

    private:

};

#endif // SAMPLEPOINT_H
