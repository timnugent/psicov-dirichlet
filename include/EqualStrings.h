// EqualStrings.h
// Kevin Karplus
// 14 Nov 1995

// EqualStrings
//	checks to see wheter two char* strings are identical.
//	Optional third argument set to one ignores case.

#ifndef EqualStrings_H
#define EqualStrings_H

#include <ctype.h>

// String comparison, with possibility of ignoring case
// Returns 1 if strings are equal, 0 otherwise
//
inline int EqualStrings(const char*a, const char*b, int IgnoreCase=0)
{   
    const char IgnoreCaseMask = ~ ('a' ^ 'A');
    for (	; *a; a++, b++)
    {    if (*a==*b) continue;		// matches
    	if (!IgnoreCase) return 0;	// mis-match
	if (((*a) ^ (*b)) & IgnoreCaseMask) 
		return 0;	// mis-match, not just case
	
	// match in all bits except the case-sensitive bit
	// do more careful check
	// Note: for brain-damaged reasons, islower and toupper take
	//	int arguments
	if (islower(*a))
	{   if(toupper(*a)==(*b)) continue;	// case-insensitve match
	    return 0;
	}
	if (! isupper(*a) || 
		tolower(*a)!=(*b)) return 0;
    }
    return (*a == *b);
}


//CHANGE LOG:
// 15 March 2004 Kevin Karplus
// 	Foo, can't remove old-style casts hidden in toupper and tolower.

#endif
