// Filenames.h
// 1 Dec 1995	(moved out of EstCommand.h)
// 2 May 1998 (moved out of Globals.h)
// Kevin Karplus

// Commands for opening files
//	has prefix and suffixes to add to filenames
//	can handle gzipped files

#ifndef Filenames_H
#define Filenames_H

#include "zfstream.h"
#include <string.h>

using namespace std;

class Filenames
{   
    // static (global) variables    
	
	static char * InFilePrefix;	// add this prefix to input file names
	static char * OutFilePrefix;	// add this prefix to output file names
	
	static char * InFileSuffix;	// add this suffix to input file names
	static char * OutFileSuffix;	// add this suffix to output file names
    
    public:	// static

	void clear(void);	// delete all suffixes and prefixes
	
	static void set_in_file_prefix(const char *fix)
	{   delete [] InFilePrefix;
	    InFilePrefix = new char[strlen(fix)+1];
	    strcpy(InFilePrefix, fix);
	}
        static void set_out_file_prefix(const char *fix)
	{   delete [] OutFilePrefix;
	    OutFilePrefix = new char[strlen(fix)+1];
	    strcpy(OutFilePrefix, fix);
	}
        static void set_in_file_suffix(const char *fix)
	{   delete [] InFileSuffix;
	    InFileSuffix = new char[strlen(fix)+1];
	    strcpy(InFileSuffix, fix);
	}
        static void set_out_file_suffix(const char *fix)
	{   delete [] OutFileSuffix;
	    OutFileSuffix = new char[strlen(fix)+1];
	    strcpy(OutFileSuffix, fix);
	}
	
	static char* full_in_filename(const char*filename);
	static char* full_out_filename(const char*filename);
	// create a newly allocated char*, 
	//	by adding the prefix and the suffix to the filename

	// Add the prefix and the suffix to the filename, 
	// and open the file for input or output.
	// If prefix fails, try without prefix.
	static gzifstream* open_input(const char*filename,
		const char* err_prefix="Error: ");
	static gzFile open_gzFile_input(const char* filename,
		const char* err_prefix="Error: ");
	static ofstream* open_output(const char*filename,
		const char* err_prefix="Error: ");

	// Open the file for input or output. (no prefix or suffix added)
	static gzifstream* open_input_fullname(const char*fullname,
		const char* err_prefix="Error: ");
	static gzFile open_gzFile_input_fullname(const char* filename,
		const char* err_prefix="Error: ");
	static ofstream* open_output_fullname(const char*filename,
		const char* err_prefix="Error: ");

	
};


// CHANGE LOG:
// 2 May 1998 Kevin Karplus
// 21 July 1997 Kevin Karplus
//	moved out of Globals
#endif
