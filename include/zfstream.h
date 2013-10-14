// zfstream.cc
// Originally by Kevin Rudland (kevin@rodin.wustl.edu)
//
// 2 Jun 1998	Kevin Karplus
//	Eliminated gzomanip, since it wasn't working in g++ to change compression
//	level, and was breaking cxx.  This eliminated any templates, and
//	templates are a major headache in c++.
//	WARNING: the compression level is hardcoded on output by Rudland,
//	so the "setcompressionlevel" function is useless.
//
// 15 August 2003 George Shackelford
//	replaced with gzstream  - simple shell class left

#ifndef _zfstream_h
#define _zfstream_h

#include <fstream>
#include "gzstream.h"
#include "zlib.h"

using namespace std;

class gzifstream : public igzstream {
public:
	gzifstream(const char* name) : igzstream(name,std::ios::in) {}
};

class gzofstream : public ogzstream {
};

/*
inline gzofstream &setcompressionlevel( gzofstream &s, int l ) 
{
  (s.rdbuf())->setcompressionlevel(l);
  return s;
}

inline gzofstream &setcompressionstrategy( gzofstream &s, int l ) 
{
  (s.rdbuf())->setcompressionstrategy(l);
  return s;
}
*/

#endif
