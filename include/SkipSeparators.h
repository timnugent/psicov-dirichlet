// SkipSeparators.cc
// Kevin Karplus
//	originally created for reading MSSP files,
//	then generalized for reading any separator-based list

#include "Input.h"
#include <assert.h>

//  Skip dest_sep-at_sep  Separators
//  return 1 if ok, 0 if not that many found.
//	reset at_sep to at_sep + number of skipped separators
//
// Normal use:  at_sep is set to 1 at the beginning of a line
//    SkipSeparators(in, at_sep, 15, ',')
//	is used to skip to the 15th field of a comma-separated list
//	(assuming you are in an earlier field).
//
//

int SkipSeparators(istream& in, int& at_sep, int dest_sep, char Separator)
{
    assert(at_sep <= dest_sep);
    in.clear();	// restore good state of input
    for (	; at_sep< dest_sep; at_sep++)
    {   char c=0;
    	while(in.get(c) && c!=Separator) {}
	if ( c!= Separator)
	    return 0;
    }
    return 1;
}


// If you don't need to keep track of where you are in the line,
// this call may be a bit simpler, since you only need to specify how many
// separators to skip.
// For example:
//    SkipSeparators(in, 1, '\n') 
//	skips to the beginning of the next line

int SkipSeparators(istream &in, int num_to_skip, char Separator)
{
    assert (num_to_skip >=0 );
    in.clear();
    for (int num=0; num< num_to_skip; num++)
    {   char c=0;
    	while(in.get(c) && c!=Separator) {}
	if ( c!= Separator)
	    return 0;
    }
    return 1;
    
}

//CHANGE LOG
// 20 Dec 1995 Kevin Karplus
// Fixed uninitialized variable in SkipSeparators
