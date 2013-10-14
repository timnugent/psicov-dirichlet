// NamedObject.cc
// 4 Dec 1995
// Kevin Karplus

#include "NamedObject.h"

#include <string.h>

const char* NamedObject::DefaultHelpString="No help available.";

void NamedObject::set_name(const char* x)
{   delete [] Name;
    Name=0;
    if (!x)	return;
    int len=strlen(x);
    Name = new char[len+1];
    strcpy(Name, x);
}

void NamedObject::read_name(istream &in)
{   char word[300];
    get_word(in,word,'=');
    if (!EqualStrings(word,"Name",1))
    {    cerr << "Error: expecting \"Name =\", but saw "
		<< word << "\n" << flush;
	return;
    }
    get_word(in, word);
    set_name(word);
}

void DeleteNamedObject(NamedObject *x)
{    delete x;
}
