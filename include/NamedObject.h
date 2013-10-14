// NamedObject.h
// Kevin Karplus
// 14 Nov 1995

// NamedObject
//	is intended to be a base class 
//	(perhaps with multiple inheritance and other base classes)
//	for objects that have names.
//
//	No symbol table is provided, though NameToPtr may
//	be appropriate to use.

#ifndef NamedObject_H
#define NamedObject_H

#include <iostream>
#include "Input.h"
#include "EqualStrings.h"


class NamedObject
{
    private:
	char *Name;
	const char *HelpString;	// NOT OWNED BY OBJECT (usually static)
	static const char* DefaultHelpString;
    
    public:
	NamedObject(void) {Name=0; HelpString=0;}
	NamedObject(const char*nm, const char*h=0)
	{   Name=0; set_name(nm);
	    set_help(h);
	}
	virtual ~NamedObject(void) {delete [] Name;}
	
	const char* name(void) const {return Name;}
	const char* help(void) const 
	{    return HelpString? HelpString: DefaultHelpString;
	}
	
	virtual void set_name(const char* x);
	virtual void set_help(const char*h) 	{    HelpString=h;	}
	
	void read_name(istream &in);
	void write_name(ostream &out) const
	{   out << "Name = " << Name << "\n";
	}
};

void DeleteNamedObject(NamedObject *x);
	// DeleteNamedOjbect just calls delete---the
	// function exists for use in NameToPtr::ApplyAll 


// A NamedToken is just a NamedOject with an extra int field for
// holding a value.

class NamedToken: public NamedObject
{
    public: 
    	int TokenValue;
	NamedToken(void): NamedObject() {};
	
	NamedToken(const char*nm, int value, const char*h=0)
	    : NamedObject(nm,h)
	{    TokenValue = value;
	}
	
};

#endif

// CHANGE LOG:
// 4 Dec 1995 Kevin Karplus
//	added help()
//	Moved set_name and read_name to new NamedObject.cc file
// 22 Dec 1995 Kevin karplus
//	added constructor from name and help strings
// 5 May 1998 Spencer Tu
//      Changed destructor to virtual so that it could be
//      deleted in a NameToPtr table with ApplyAll
// 14 May 2002 Kevin Karplus
//	Added NamedToken type.
//	set_name() made virtual, so that derived objects
//	can redefine it.
