// File: NameToPtr.cc
// Programmer: Kevin Karplus
// Creation Date: 11 Jan 1994

// Abstract:
//	NameToPtr is a hash table for mapping the name of a thing to a
//		pointer to the thing.
//	For versatility, it returns a void * pointer---and so the NameToPtr
//		class objects should always be buried as members of some
//		class that knows what sorts of objects it plans to keep
//		in the table.

#include "NameToPtr.h"
#include "EqualStrings.h"
#include <assert.h>
#include <string.h>	// for strlen
#include <iostream>	// for cerr
#include <stdlib.h>	// for abort()


NameToPtr::NameToPtr(int size)
{   assert(size>0);
    HashSize = size;
    NameTable = new NameToPtrItem* [size];
    NumNames=0;
    for (int i=size-1; i>=0; i--)
	NameTable[i] = 0;
    ignore_case(0);
}


NameToPtr::~NameToPtr()
{
    for (int i=HashSize-1; i>=0; i--)
    {	if (NameTable[i]) 
    	    delete NameTable[i];
    }
    delete [] NameTable;
}

// This is a crude but fast hash function.
// There are several better ones that could be used if this one turns
// out to be inadequate.
//
const char IgnoreCaseMask = ~ ('a' ^ 'A');

int NameToPtr::Hash(const char *n) const
{
    unsigned int h=0;
    for(	;*n; n++)
        h = (h<<3) ^ ((*n) & IgnoreCaseMask);
    return h % HashSize;	// must be in range 0..HashSize-1
}


void NameToPtr::Rehash(int NewSize)
{   assert(NewSize>0);
    NameToPtrItem **oldtable=NameTable;
    int oldsize=HashSize;
    
    NameTable = new NameToPtrItem*[NewSize];
    HashSize = NewSize;
    for (int i=oldsize-1; i>=0; i--)
    {	NameToPtrItem* nextn=0;
        for (NameToPtrItem* n= oldtable[i]; n; n=nextn)
	{   int h=Hash(n->object->name());
	    nextn=n->next;
	    n->next=NameTable[h];
	    NameTable[h] = n;
	}
    }
    delete [] oldtable;
}

NamedObject* NameToPtr::FindOldName(const char * name, OptionIfNew ifnew) const
{
    assert (ifnew != CreateIfNew);
    int h = Hash(name);
    for (NameToPtrItem *n=NameTable[h]; n; n=n->next)
    {   if (EqualStrings(n->object->name(), name, IgnoreCase))
    		return n->object;
    }
    switch (ifnew)
    {case ErrorIfNew:
    	cerr << "NameToPtr::FindOldName expected to find name " << name 
		<< " in the table of " << NumNames << " names.\n" << flush;
	abort();
    case CreateIfNew:  // shouldn't happen but the compiler wants it
	abort();
	break;

    case ZeroIfNew:
    	return 0;
    }
    assert(0);	// unreachable code
    return 0;
}

void NameToPtr::AddName(NamedObject* object, 
	OptionIfOld ifold)
{
    const char* name = object->name();
    int h = Hash(name);
    NameToPtrItem *n;
    for (n=NameTable[h]; n; n=n->next)
    {   if (EqualStrings(n->object->name(), name, IgnoreCase))
    	{   if (ifold==ErrorIfOld)
	    {	cerr << "Error: NameToPtr::AddName "
			<< "already has an object named " 
			<< n->object->name() << " so can't add "
			<< name
			<< "\n" << flush;
		abort();
	    }
	    return;
	}
    }
    n = new NameToPtrItem(object);
    n->next = NameTable[h];
    NameTable[h] = n;
    NumNames ++;
}

void NameToPtr::DeleteName(const char *name, OptionIfNew ifnew)
{
    if (name==NULL) return;
    assert(ifnew!=CreateIfNew);		// zero or error only
    int h = Hash(name);
    NameToPtrItem* oldn=0;
    for (NameToPtrItem *n=NameTable[h]; n; oldn=n, n=n->next)
    {   if (EqualStrings(n->object->name(), name, IgnoreCase))
    	{   if (oldn) oldn->next=n->next;
	    else NameTable[h] = n->next;
	    n->next=0;
	    delete n;
	    NumNames--;
	    return;
	}
    }
    if (ifnew==ErrorIfNew)
    {	cerr << "NameToPtr::DeleteName expected to find name " << name
		<< " in the table of " << NumNames << " names.\n" << flush;
	abort();
    }
}


void NameToPtr::apply_all(void (* fun)(NamedObject *))
{
    for(int h=HashSize-1; h>=0; h--)
    {   for (NameToPtrItem *n=NameTable[h]; n; n=n->next)
    	{    fun(n->object);
	}
    }
}

void NameToPtr::apply_all(void (* fun)(const NamedObject *)) const
{
    if (NumNames==0) return;
    for(int h=HashSize-1; h>=0; h--)
    {   for (NameToPtrItem *n=NameTable[h]; n; n=n->next)
    	{    fun(n->object);
	}
    }
}


// CHANGE LOG:
// 13 January 1994, Kevin Karplus
//	Initialized hashtable to all zeros.
// 13 Nov 1995 Kevin Karplus
//	Fixed bug in EqualStrings (wasn't checking case bit right)
// 14 Nov 1995 Kevin Karplus
//	Moved EqualStrings to separate .h file
//	Modified code throughout to use NamedObject
//	Now keep names entirely in the objects.
// 8 August 2003 Kevin Karplus
//	Added NumNames==0 check to apply_all, to speed up empty applies.
// 21 July 2004 Kevin Karplus
//	Added check to avoid trying to delete empty name pointers.
