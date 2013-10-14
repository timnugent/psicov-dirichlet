// log2.h	
// Kevin Karplus
// 31 May 1995
//	Should be in math.h

#ifndef LOG2_H
#define  LOG2_H


#include <math.h>


//inline double log2(double x) 
// {    return M_LOG2E *log(x);
//}

inline double clip_log2(double p) {return p<=0.? -1000.: log2(p); }

// Thu Jun 30 12:00:34 PDT 2005 Kevin Karplus
// Commented out inline log2, since it seems to be in modern
//	g++ math libraries
//	(note this may need checking on different systems)
#endif
