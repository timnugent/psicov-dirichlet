#ifndef abs_H
#define abs_H

// abs function now in cmath - 17 July 2003 G. Shackelford

// Thu Jan 26 15:23:11 PST 2006 Kevin Karplus
//	g++ 3.4.4 doesn't handle abs well, so this is needed again.

namespace kk {

// get int abs from stdlib with C syntax?

inline short abs(short x)	{return x<0? -x : x;}
inline int abs(int x)		{return x<0? -x : x;}
inline long abs(long x)		{return x<0? -x : x;}
inline float abs(float x)	{return x<0? -x : x;}
inline double abs(double x)	{return x<0? -x : x;}

}



// 15 March 2004 Kevin Karplus
//	moved inline to beginning of line.

#endif
