// Filenames.cc
// Kevin Karplus
// 22 Nov 1995

#include <iostream>
#include <fstream>
#include <stdio.h>	// for popen
#include <string.h>	// for strcat
#include <stdlib.h>	// for system
#include <assert.h>

#include "Filenames.h"

char * Filenames::InFilePrefix=0;
char * Filenames::OutFilePrefix=0;
char * Filenames::InFileSuffix=0;
char * Filenames::OutFileSuffix=0;

void Filenames::clear(void)
{
    delete [] InFilePrefix;    InFilePrefix=0;
    delete [] OutFilePrefix;   OutFilePrefix=0;
    delete [] InFileSuffix;    InFileSuffix=0;
    delete [] OutFileSuffix;    OutFileSuffix=0;
}

// create a newly allocated char* string with the full filename
//	(adding prefix and suffix, if any)
char * Filenames::full_in_filename(const char* filename)
{
    assert(filename!=0);
    int len = strlen(filename);
    int pre_len= InFilePrefix? strlen(InFilePrefix): 0;
    int suf_len= InFileSuffix? strlen(InFileSuffix): 0;
    
    char* full_name= new char[pre_len+len+suf_len+1];
    full_name[0]=0;
    if (pre_len) strcpy(full_name, InFilePrefix);
    strcat(full_name, filename);
    if (suf_len) strcat(full_name, InFileSuffix);
    return full_name;
}
// create a newly allocated char* string with the full filename
//	(adding prefix and suffix, if any)
char * Filenames::full_out_filename(const char* filename)
{
    assert(filename!=0);
    int len = strlen(filename);
    int pre_len= OutFilePrefix? strlen(OutFilePrefix): 0;
    int suf_len= OutFileSuffix? strlen(OutFileSuffix): 0;
    
    char* full_name= new char[pre_len+len+suf_len+1];
    full_name[0]=0;
    if (pre_len) strcpy(full_name, OutFilePrefix);
    strcat(full_name, filename);
    if (suf_len) strcat(full_name, OutFileSuffix);
    return full_name;
}

gzifstream* Filenames::open_input_fullname(const char*full_name,
	const char* err_prefix)
{
    gzifstream* file= new gzifstream(full_name);
    if (!file || !file->good())
    {   char* gzip_name= new char[strlen(full_name) +4];
        strcpy(gzip_name, full_name);
	strcat(gzip_name, ".gz");
	delete file;
	file = new gzifstream(gzip_name);
	if (!file || !file->good())
        {   if (err_prefix!=NULL)    cerr << err_prefix;
	    cerr << "Couldn't open file " << full_name 
		<< " or " << gzip_name
		<< " for input\n";
            delete file;
	    file=0;
	}
	delete [] gzip_name;
    }
    return file;
}

gzifstream* Filenames::open_input(const char*filename,
	const char* err_prefix)
{   char *full_name=full_in_filename(filename);    
    gzifstream *file = open_input_fullname(full_name, "Warning: ");
    delete [] full_name;
    if (!file)
    {	cerr << "Trying " << filename << "\n"
    		<< flush;
        file = open_input_fullname(filename, err_prefix);
    }
    return file;
}



gzFile Filenames::open_gzFile_input_fullname(const char*full_name,
	const char* err_prefix)
{
    gzFile file= gzopen(full_name,"r");
    if (!file)
    {   char* gzip_name= new char[strlen(full_name) +4];
        strcpy(gzip_name, full_name);
	strcat(gzip_name, ".gz");
	file = gzopen(gzip_name, "r");
	if (!file)
        {   if (err_prefix!=NULL)    cerr << err_prefix;
	    cerr << "Couldn't open file " << full_name 
		<< " or " << gzip_name
		<< " for input\n";
	}
	delete [] gzip_name;
    }
    return file;
}

gzFile Filenames::open_gzFile_input(const char*filename,
	const char* err_prefix)
{   char *full_name=full_in_filename(filename);    
    gzFile file = open_gzFile_input_fullname(full_name, "Warning: ");
    delete [] full_name;
    if (!file)
    {	cerr << "Trying " << filename << "\n"
    		<< flush;
        file = open_gzFile_input_fullname(filename, err_prefix);
    }
    return file;
}



ofstream* Filenames::open_output_fullname(const char*full_name,
	const char* err_prefix)
{
    ofstream* file= new ofstream(full_name);
    if (!file || !file->good())
    {   if (err_prefix!=NULL)    cerr << err_prefix;
        cerr << "Couldn't open file " << full_name << " for output\n";
	delete file;
	return 0;
    }
    return file;
}

ofstream* Filenames::open_output(const char*filename,
	const char* err_prefix)
{   char *full_name=full_out_filename(filename);    
    ofstream *file = open_output_fullname(full_name, "Warning: ");
    delete [] full_name;
    if (!file)
    {    file = open_output_fullname(filename, err_prefix);
    }
    return file;
}



// CHANGE LOG:
// 21 July 1997 Kevin Karplus
//	copied from estimat-dist and removed everything related to
//	TestSet or regularizer stack.
// 23 July 1997 Kevin Karplus
//	modified open_input to try gunzipping file if the file
//	doesn't exist but the file.gz does
// 14 Sept 1997 Kevin Karplus
//	FINALLY fixed open_input to use popen and zcat if need to gunzip
// 22 March 1998 Kevin Karplus
//	Repaired damage Christian did to open_input, by providing
//	option of adding prefix and suffix or not.
// 22 May 1998 Kevin Karplus
//	Moved out of Globals
// 3 June 1998 Kevin Karplus
//	Modified to use zfstream.
// 29 March 2002 Kevin Karplus
//	Modified open_input and open_output to try without
//	prefix if open fails with prefix.
// 30 Aug 2002 Kevin Karplus
//	Added err_prefix arguments to all open routines,
//	so that error messages are easier to find in output.
// 17 July 2003 George Shackelford
//	Added and changed include files for gcc 3.x
// 23 June 2005 Kevin Karplus
//	Took out extra printout of err_prefix in open_input before
//	"Trying" message
