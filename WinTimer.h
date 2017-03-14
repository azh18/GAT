#pragma once
#include<windows.h>

class MyTimer
{
private:
	LARGE_INTEGER _freq;
	LARGE_INTEGER _start;
	LARGE_INTEGER _stop;
public:

	MyTimer()
	{
		QueryPerformanceFrequency(&_freq);
	}

	inline void start()
	{
		QueryPerformanceCounter(&_start);
	}

	inline void stop()
	{
		QueryPerformanceCounter(&_stop);
	}

	inline float elapse()
	{
		return 1e3*(_stop.QuadPart - _start.QuadPart) / _freq.QuadPart;
	}

	inline long long ticks()
	{
		return _stop.QuadPart - _start.QuadPart;
	}
};