// File: NameToPtr.h
// Programmer: Kevin Karplus
// Creation Date: 11 Jan 1994

// Abstract:
//	NameToPtr is a hash table for mapping the name of a thing to a
//		pointer to the thing.
//	The thing must be derived from the NamedObject class.
//	Since the NamedObject knows nothing about the hash table,
//	it is an error to change the name of an object while it is
//	in the table.
//

// To Do:
//	Should the OptionIfNew and OptionIfOld enums be moved into
//	NameToPtr, so that one has to use NameToPtr::ZeroIfNew?

#ifndef NAMETOPTR_H

#include "NamedObject.h"

// A NameToPtrItem is an element of a linked list of NamedObjects.
//	The object is not owned by list element.
//
class NameToPtrItem
{   friend class NameToPtr;
    
    private:
	NamedObject *object;
		// although NameToPtr never changes an object,
		// it needs to be able to return objects, not
		// just const objects, so can't use const.
	NameToPtrItem *next;
    
	NameToPtrItem(NamedObject *o)
	{    object=o;     next=0;	    }
    
	~NameToPtrItem()	{delete next;}
};


// There are three standard ways to handle not finding a name in a table:
//	ZeroIfNew	return 0, so that the calling program knows
//			the name is new. 
//	ErrorIfNew	report an error and abort---the calling
//			program expects the name to be there.
//	CreateIfNew	If the name is not already there, create a new
//			object with the name and put it there.
//			The "CreateIfNew" option just adds the name to the
//			hashtable with a zero pointer, since it doesn't
//			know what sort of object to create.
//
	typedef enum {ZeroIfNew, ErrorIfNew, CreateIfNew} OptionIfNew;

// When looking up a name, we may also want to be sure that the name
//	does not previously exist--there are two options here
//	OKIfOld		just return the old thing
//	ErrorIfOld	shouldn't already be something with this name.
//	ReplaceIfOld	simply replace the old pointer with a new one
//			WARNING: old object not deleted (don't know type)
//			This could result in massive storage leakage if
//			NameToPtr owns the objects.
	typedef enum {OKIfOld, ErrorIfOld, ReplaceIfOld} OptionIfOld;


// case may be considered or ignored on any operation.
// Names are stored with original case, and hashing function ignores case,
//	so the case sensitivity can be turned on and off at will.

class NameToPtr
{   
    	int HashSize;	// number of initial slots in hashtable.
    	int NumNames;	// number of names in the hashtable
	NameToPtrItem **NameTable;	// chained-overflow hashtable
	
	int IgnoreCase;		// default is to consider case relevant
	
	int Hash(const char *name) const; // map the string to [0..HashSize-1]
	
    public:
	NameToPtr(int size=937);
	void Rehash(int NewSize);	// grow (or shrink) the hashtable
	// for best hashing, size should be a prime number
	
	~NameToPtr();
	
	void ignore_case(int ignore=1)
	{   IgnoreCase=ignore;
	}
	
	int RetNumNames() const {return NumNames;}
	int RetHashSize() const {return HashSize;}
	
	NamedObject* FindOldName(const char * name, 
			OptionIfNew ifnew=ErrorIfNew) const;
        // Note: CreateIfNew is not a legal option---it is provided
	//	for use by higher level programs that may know
	//	how to create objects.
	
	void AddName(NamedObject *object, 
			OptionIfOld ifold=ReplaceIfOld);
	//	Since objects may be searched for with or without case,
	//	it is an error to enter two names which differ only in
	//	the case of the names when ignore_case is set.

	
	void DeleteName(const char *name, OptionIfNew ifnew=ErrorIfNew);
	// ZeroIfNew causes a no-op if the name is not found.
	
	void apply_all(void (*fun)(NamedObject*));	// apply function to all 
	void apply_all(void (*fun)(const NamedObject*)) const;	// apply function to all 
};

// Deleting a NameToPtr table does not delete the NamedObjects the
// table includes.
// If you want to delete them, start by calling ApplyAll(&DeleteNamedObject)

#define NAMETOPTR_H (1)
#endif

// Change Log:
// 14 Nov 1995 Kevin Karplus
//	Changed to return NamedObject* instead of void*
//	Also, names now stored exclusively in the objects, not copied
//	into the hashtable.
//	Eliminated FindName, modified ApplyAll to apply to NamedObjects.
//	Modified AddName to get name from object.
// 26 Nov 1996 Spencer Tu
//      Changed definition of NameToPtr::FunctionNameObj
//      from
//
//          typedef void FunctionNameObj (NamedObject*);
//
//      to
//
//          typedef void (* FunctionNameObj)(NamedObject*);
//
//      so cxx would accept passing an argument of the type
//      to a function. Same for NameToPtr::FunctionNameConstObj .
// 24 Nov 1999
//	Renamed ApplyAll to apply_all and eliminated FunctionNameObj type
