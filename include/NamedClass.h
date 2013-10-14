// NamedClass.h
// Kevin Karplus
// 13 Nov 1995

// NamedClass
//	is intended for use as base class for classes that need to
//	have run-time type determination.
//
// IdObject
//	is the type used for giving unique IDs to classes, and
//	and for maintaining an "isa" heirarchy.  For ease of use,
//	the name lookup ignores case.
//
//	The class name "0" can be used for generating null pointers.

#ifndef NAMEDCLASS_H
#define NAMEDCLASS_H

#include "NamedObject.h"
#include "NameToPtr.h"
#include <iostream>
#include <vector>


class NamedClass;	// forward declaration


class IdObject: public NamedObject
{
    public:
         typedef NamedClass * (*NamedClassFunction) (void);
    
    private:
	NamedClassFunction Create;	// a function that creates
				// a new object of the class
				// may be 0 for pure virtual classes
	
	typedef void (*VoidFunction)(IdObject *);
	VoidFunction InitIsA;	// initialization function for IsA
				// called first time IsA is called
	
	
	int NumIsA, AllocIsA;
	 IdObject ** IsAArray;	// array of classes that this
		// class is a specialization of (transitive relation)
	
	int NumIsNotA, AllocIsNotA;
	 IdObject ** IsNotAArray;	// array of classes that this
					// class is known NOT to be 
	
	void Realloc(int NewAllocIsA, int NewAllocIsNotA);

	static NameToPtr* LookupTable;	
		// symbol table for name-to-IdObject mapping
    public:

	typedef void (*IdObjectFunction)(IdObject *i);
	
	IdObject(const char *s,			//name
		NamedClassFunction c=0,		// create function
		IdObjectFunction init=0,	// IsA initialization
		const char* hlp=0		// help string
		);		

	// Note: any pointers to this in the is_a caches of other
	// IdObjects are NOT cleaned up.  
	// IdObjects are not designed for dynamic deletion---this
	// code is only here for cleaning up memory at the end of
	// an execution.
	~IdObject(void)
	{   delete [] IsAArray;
	    delete [] IsNotAArray;
	    LookupTable->DeleteName(name(), ErrorIfNew);
	}
	
	static IdObject *id(const char * name)
	{    return dynamic_cast<IdObject*>(LookupTable->FindOldName(name, ZeroIfNew));
	}
	
	static const NameToPtr* lookup_table(void) {return LookupTable;}
	
	NamedClass * create(void) const {return Create? (*Create)(): 0;}
	
	bool is_a(IdObject *x);	// return 1 if this is an x
	void add_is_a(IdObject *x);	// record that this is an x
	
	static void apply_all(void (* fnc)(IdObject *))
	{    LookupTable->apply_all(reinterpret_cast<void (*)(NamedObject *)>(fnc));
	}
	static void apply_all(void (* fnc)(NamedObject *))
	{    LookupTable->apply_all(fnc);
	}
};

typedef IdObject* IdObjectPtr;

class NamedClass
{
    // members that have to be provided in all derived classes
    // overloading the ones in the base class:
    private:
	static IdObject ID;
	static void init_is_a(IdObject*) {}	// may be omitted if no-op

    protected:
	virtual int read_knowing_type(istream&in) {return 0;}
	// returns 1 if it has already read the terminal
	//	EndClassName = command
	
        virtual void write_knowing_type(ostream &out) const {}
	
    public:
	static  IdObject *classID(void) {return &NamedClass::ID;}
	virtual IdObject *type(void) const {return &NamedClass::ID;}
	// IdObject not const, since is_a() needs to cache results

	NamedClass(void) {}	// pure virtuals don't need a constructor
				// all others must have a void constructor
        virtual ~NamedClass() {}    
	
    // shared members for all NamedClass derived classes:
    public:
	int is_a(IdObject *x) const	{    return type()->is_a(x);  }
	
	ostream& write(ostream &out) const
	{   out << "ClassName = " << type()->name() << "\n";
	    write_knowing_type(out);
	    out << "EndClassName = " << type()->name() << "\n";
	    return out;
	}
	
	istream& read(istream &in);
	
	static NamedClass* read_new(istream &in);
		// read an object from istream, and return it (newly allocated)
		// Ojbect must start with "CLASSNAME=" (case irrelevant)
	static NamedClass* read_new_after_classname(istream &in);
		// read an object from istream, and return it (newly allocated)
		// Assumes that CLASSNAME= has already been read
    	
	static int read_many_new(vector<NamedClass*>& objects,
		istream &in,
		IdObject* OnlyThisType=&NamedClass::ID);
		// Read many objects from istream (until EOF)
		//	and add newly allocated objects to the end of
		//	the vector
		// Return the number added.
		// If OnlyThisType is set, an objects that is not of the
		//	specified type will cause and error.

	static int read_many_new(vector<NamedClass*>& objects,
		const char* filename,
		IdObject* OnlyThisType=&NamedClass::ID);
		// opens filename for input and calls
		// read_many_new(objects, istream&, OnlyThisType)
};

inline ostream& operator << (ostream&out, const NamedClass& n)
{    return n.write(out);
}

// Note: this operator requires an exact type match---use read_new() 
// if you want type to be read from input stream.
inline istream& operator >> (istream&in, NamedClass &n)
{    return n.read(in);
}

// A class for defining InputCommands for reading NamedClass ojbects
class NamedClassCommand: public NamedObject
{    private:    
       typedef int (*ICF) (istream &in, NamedClass *change, NamedClassCommand *self);
       ICF CommandFunction;
    
    public:
        NamedClassCommand(const char *nm, ICF icf)
        {   set_name(nm);
            CommandFunction = icf;
        }

        virtual int execute(istream &in, NamedClass *change)
        {   return (*CommandFunction) (in, change, this);
        }
};


// Useful generic commands:

// skip to end of line
int ReadComment(istream &in,
		    NamedClass *change,
		    NamedClassCommand *self);

// if change is also a NamedObject, then set its name.
int ReadName(istream &in,
		    NamedClass *change,
		    NamedClassCommand *self);


// Called for both the ClassName and EndClassName keywords
int VerifyClassName(istream &in, 
	    NamedClass *change,
	    NamedClassCommand* self);



// USE OF NamedClass

// Here is an example of a declaration of a pure virtual type and
// a real class derived from it.  The virtual type "virt" has no init_is_a
// function, since it is not derived from any higher class type.
// It also has no constructor, since it is pure virtual, though it could
// have one for use by derived classes, if needed.
//

// pure virtual class in .h file:
//
//  class virt: public NamedClass
//  {
//	private:
//		static IdObject ID;
//		virtual int read_knowing_type(istream &in);
//		virtual void write_knowing_type(ostream &out) const;
//	public:
//		static  IdObject *classID(void) {return &virt::ID;}
//		virtual IdObject *type(void) const {return &virt::ID;}
//		
//	// ... and whatever is needed for the class itself.
//  };

// pure virtual classes in .cc file:
//
//    IdObject virt::ID ("virt");
// OR
//    IdObject virt::ID ("virt",0,0,"a virtual class used for demo purposes");



// derived class in .h file
//
//
//  class deriv: public vert
//  {
//	private:
//		static IdObject ID;
//		static void init_is_a(IdObject*self)
//		{    self->add_is_a(virt::classID());
//		}
//		virtual int read_knowing_type(istream &in);
//		virtual void write_knowing_type(ostream &out) const;
//	public:
//		static  IdObject *classID(void) {return &deriv::ID;}
//		virtual IdObject *type(void) const {return &deriv::ID;}
//		deriv(void) {...}
//		
//	// ... and whatever is needed for the class itself.
//  };

// derived class in .cc file:
//
//    static NamedClass *create_deriv(void) {return new deriv;}
//    IdObject deriv::ID ("deriv", create_deriv, deriv::init_is_a,
//	"a derived class of virt, used for demo purposes");

// CHANGE LOG:
// 3 Dec 1995 Kevin Karplus
//	Added apply_all().
// 4 Dec 1995 Kevin Karplus
//	added help string to constructor
// 12 Jan 1996 Kevin Karplus
//	minor change to IdObject constructor declaration, to eliminate public
//	use of VoidFunction.  No change to usage.
// 12 Oct 1996 Spencer Tu
//      change def of IdObject:IdObjectFunction from
//
//          typedef void IdObjectFunction (IdObject *i);
//
//      to
//     
//          typedef void (*IdObjectFunction)(IdObject *i);
//
//      so cxx would accept it on compiling the initialization
//      of deriv::ID.
// 26 Nov 1996 Kevin Karplus
//	Fixed corresponding uses of IdObjectFunction to remove extra *
// 26 Nov 1996 Spencer Tu
//      Changed uses of NameToPtr::FunctionNameObj to remove extra *.
//      Similar to the above problem (see NameToPtr.h).
//      Affected function was IdObject::apply_all().
// 21 May 1998 Spencer Tu
//      Added a virtual destructor for NamedClass so that NamedClass
//      objects can be deleted when they are put in a NameToPtr table
//      with a call to NameToPtr::ApplyAll()
// 5 Nov 1999 Kevin Karplus
//	split read_new_after_classname off from read_new
// 27 Feb 2004 Kevin Karplus
//	Converted old-style casts.
// Sat Jun 18 09:56:02 PDT 2005 Kevin Karplus
//	Added read_many_new() static members.
// Wed Oct 25 10:27:31 PDT 2006 Kevin Karplus
//	Made read_knowing_type and write_knowing_type protected,
//	instead of private. 
// Wed Oct 25 10:29:49 PDT 2006 Kevin Karplus
//	Added read() for already allocated objects (insists on typematch).
//	Made write() return ostream&.
// Mon Dec 10 16:53:13 PST 2007 Kevin Karplus
//	Removed extraneous NamedClass:: inside declaration of NamedClass
// Sun Apr 19 17:05:55 PDT 2009 Kevin Karplus
//	Added NamedClassCommand type.
#endif
