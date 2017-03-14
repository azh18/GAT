#include "MBB.h"



MBB::MBB()
{
	xmin = 0;
	xmax = 0;
	ymin = 0;
	ymax = 0;
}

MBB::MBB(float val_xmin,float val_ymin,float val_xmax,float val_ymax)
{
	xmin = val_xmin;
	xmax = val_xmax;
	ymin = val_ymin;
	ymax = val_ymax;
}

bool MBB::pInBox(float x, float y) {
	if (x <= xmax &&x >= xmin&&y <= ymax&&y >= ymin)
		return true;
	else
		return false;
}

int BoxIntersect(MBB& a1, MBB& b1) {
	MBB a, b;
	bool swaped = false;
	if (a1.xmin < b1.xmin) {
		a = a1;
		b = b1;
	}
	else if (a1.xmin == b1.xmin) {
		if (a1.xmax > b1.xmax) {
			a = a1;
			b = b1;
		}
		else
		{
			b = a1;
			a = b1;
			swaped = true;
		}
	}
	else
	{
		b = a1;
		a = b1;
		swaped = true;
	}
	if (b.xmin >= a.xmax)
		return 0;
	else
	{
		if (b.ymax <= a.ymin)
			return 0;
		else if (b.ymin >= a.ymax)
			return 0;
		else {
			if (a.pInBox(b.xmin, b.ymin) && a.pInBox(b.xmin, b.ymax) && a.pInBox(b.xmax, b.ymin) && a.pInBox(b.xmax, b.ymax))
			{
				if (!swaped)
					return 2;
				else
					return 3;
			}
			else
				return 1;
		}
	}
}
/* return 0:不相交
return 1:相交但不包含
return 2:a1包含b1
return 3:b1包含a1
*/

int MBB::intersect(MBB& b) {
	return (BoxIntersect(*this, b));
}
/* return 0:不相交
   return 1:相交但不包含
   return 2:this包含b
   return 3:b包含this
*/



MBB::~MBB()
{
}
