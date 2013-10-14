// Input.h
// Kevin Karplus
// 8 August 1994

// ABSTRACT
// This file declares various useful input routines, providing uniform
// ways for dealing with comments in input files.
//
//
//	void get_word(istream &, char* word)
//		reads a (white-space-terminated) word from the input
//		into a pre-allocated buffer provided.
//		Ignores comments starting with # or // 
//		WARNING: no length check is done---this should only
//		be used for known-safe input.
//
//	void get_word(istream &, char* word, char stop_char, int skip)
//		reads a (white-space-terminated) word from the input
//		into a pre-allocated buffer provided.
//		Ignores comments starting with # or // 
//		Stops immediately *after* reading stop_char.
//		if the word is terminated by white-space and not stop_char,
//		   then  {if skip,  SkipSeparators is called to skip through
//					the next stop_char,
//		          else skip white space (and stop_char if present)}
//		A comment is treated as a '\n' character for stopping.
//		WARNING: no length check is done---this should only
//		be used for known-safe input.
//
//	bool verify_word(istream&, const char* word)
//		reads a word and makes sure it is identical to the
//		specified word
//
//	const char *read_type(istream&)
//		reads a word and returns the canonical classID for the
//		class that has registered that word as its name.
//		If no class has registered the name, read_type returns 0.
//
//	int SkipSeparators(istream& in, 
//			int& at_sep, int dest_sep, char Separator)
//		Skip dest_sep-at_sep  Separators
//  		return 1 if ok, 0 if not that many found.
//		reset at_sep to at_sep + number of skipped separators
//
// 		Normal use:  at_sep is set to 1 at the beginning of a line
//    			SkipSeparators(in, at_sep, 15, ',')
//		is used to skip to the 15th field of a comma-separated list
//		(assuming you are in an earlier field).
//
//	int SkipSeparators(istream &in, int num_to_skip, char Separator)
//		If you don't need to keep track of where you are in the line,
// 		this call may be a bit simpler, since you only need to
//		specify how many separators to skip.
// 		For example:
//    			SkipSeparators(in, 1, '\n') 
//		skips to the beginning of the next line

#ifndef INPUT_H
#define INPUT_H

#include <iostream>

using namespace std;

void get_word(istream &in, char *word);  
char get_word(istream &in, char *word, char stop_char, bool skip=1);
bool verify_word(istream &in, const char *should_be);


char *read_line(istream&in);
// Read a line (skipping comments and blank lines)
//	and strip leading white space.
// Allocate new storage for it, and return the char* string.
// Return 0 if no more data.

int SkipSeparators(istream& in, int& at_sep, int dest_sep, char Separator);
int SkipSeparators(istream &in, int num_to_skip, char Separator);

//CHANGE LOG:
// Fri Feb 11 06:31:15 PST 2005 Kevin Karplus
//	Made skip parameter of get_word be bool.
#endif
