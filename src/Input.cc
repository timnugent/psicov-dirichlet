// Input.cc
// Kevin Karplus
// 8 August 1994


// TO DO:
//	Redo verify_word so that it doesn't need get_word (check a
//	character at a time).

#include "Input.h"
#include <string.h>

// read a word, ignoring white space 
//	comments starting with //	end with \n
//	comments starting with #	end with \n
//
void get_word(istream &in, char *word)
{   
    if (! in.good())
    {   word[0]=0;
    	return;
    }
    for(in>>word;
	in.good() && (word[0]=='/' && word[1]=='/' || word[0]=='#') ;
	in>>word)
    {
        char c;	
        while( in.get(c) && c!='\n')
        {}
    } 
    if (!in.good())
	word[0] = 0;
}


// similar, but with extra option for early termination,
//	even if no word found
// read a word, ignoring white space 
//	comments starting with //	end with \n
//	comments starting with #	end with \n
//
// If stop_char is non_zero, stop after stop_char is read (as well as stopping
//	for white space)
// If skip is set, then skip over everything until stop_char is read.
//
// Usual use: stop_char is '\n' to avoid reading past end of line.
// Note: comments are treated as '\n' 
//
// return last character read
//
char get_word(istream &in, char *word, char stop_char, bool skip)
{   
    word[0] = 0;
    if (! in.good())   return 0;
    char *w=word;    
    char c=0;
    while (word[0] == 0)	// haven't read anything yet
    {    // skip white-space
	for (in.get(c);
	    in.good() && c!=stop_char 
	    	&& (c==' ' || c=='\t' || c=='\n' || c=='\f');
	    in.get(c))
	{}
	
	if (!in.good()) return c;
	if (c==stop_char) return c;

	if (c=='#') 		// skip # comment
	{   while( in.get(c) && c!='\n')
	    {}
	    if (stop_char == '\n') return c;
	    else continue;
	} 
	else if (c=='/')
	{   if (in.get(c) && c=='/')		// skip // comment
	    {   while( in.get(c) && c!='\n')
		{}
		if (stop_char == '\n') return c;
		else continue;
	    }
	    else
	    {   *w++ = '/';
	        *w = 0;
	    }
	}
	
	*w++ = c;
	*w = 0;
	while(in.get(c) && (c!=stop_char 
			&& c!=' ' && c!= '\t' && c!='\n' && c!='\f'))
	{   *w++ = c;
	    *w = 0;
	}
	if (c==stop_char && !skip)
	{    in.putback(c);
	}
	
	if (stop_char && c!=stop_char) 
	{   if (skip)
            {
                SkipSeparators(in, 1, stop_char);
                c = stop_char;
            }
	    else
	    {
                char oldc = c;
                for (;
		    in.good()&& c!=stop_char 
		    	&& (c==' ' || c=='\t' || c=='\n' || c=='\f');
		    in.get(c))
		{
                    oldc = c;
                }
                in.putback(c);
	    }
	}
    }

    return c;
}

bool verify_word(istream &in, const char *should_be)
{
    char word[100];
    if (!in.good() && should_be[0]!=0)	return 0;
    get_word(in,word);
    if (strcmp(word,should_be)==0)
    	return 1;
    cerr << "Error: expecting " << should_be << " saw " << word 
    	<< "\n" << flush;
    return 0;    
}

// CHANGE LOG:
// 1 Sept 1995 Kevin Karplus
//	Added new get_word with stop_char
// 17 Nov 1995 Kevin Karplus
//	Modified get_word with stop_char to call SkipSeparators
//	to make sure stop_char is swallowed.
// 21 Feb 1996 Kevin Karplus
//	Added "skip" variable to get_word to control whether the
//	the routine skipped until it had swallowed the stop_char or not.
// 1 June 1996 Spencer Tu
//      Modified get_word() to return last char read.
// 12 Oct 1996 Spencer Tu
//      In get_word(), replaced call to istream::unget() with a call
//      to istream::putback().  Unget is not supported by cxx lib,
//      but putback() is supported by both cxx and gnu.
// 27 December 1996 Kevin Karplus
//	Changed meaning of skip==0 in get_word, so that the stop
//		character is NOT swallowed.
// 23 August 2004 Kevin Karplus
//	Added extra in.good() tests before reading words.
//	Changed verify_word to return bool.
// Wed Jun 29 06:50:24 PDT 2005 Kevin Karplus
//	Modified get_word so that formfeeds are treated as white space,
//	and in.good checked before value of c in returns.
