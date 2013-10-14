// NamedClass.cc
// Kevin Karplus
// 13 Nov 1995

#include "NamedClass.h"
#include "Input.h"
#include "EqualStrings.h"
#include "Filenames.h"
#include <assert.h>
#include <stdlib.h> //for abort()

NameToPtr* IdObject::LookupTable=0;

// ZeroObject is provided so that "ClassName = 0" "EndClassName = 0"
// can be used to indicate empty pointers.

static IdObject ZeroObject("0");

IdObject NamedClass::ID("NamedClass");


IdObject::IdObject(const char *s, NamedClassFunction c,
	VoidFunction init,
	const char* hlp)
{   if (!LookupTable)	
    {	LookupTable = new NameToPtr(59);
	LookupTable->ignore_case(1);
    }
    set_name(s);
    set_help(hlp);
    Create = c;
    InitIsA = init;
    NumIsA=AllocIsA=NumIsNotA=AllocIsNotA=0;
    IsAArray=IsNotAArray=0;
    LookupTable->AddName(this, ErrorIfOld);
}


void IdObject::Realloc(int NewIsSize, int NewIsNotSize)
{
    if (NewIsSize != AllocIsA)
    {   assert(NewIsSize >= NumIsA);
	AllocIsA = NewIsSize;
	 IdObject** tmp = new IdObjectPtr[AllocIsA];
	for (int i=NumIsA-1; i>=0; i--)
	    tmp[i] = IsAArray[i];
	delete [] IsAArray;
	IsAArray = tmp;
    }
    if (NewIsNotSize != AllocIsNotA)
    {   assert(NewIsNotSize >= NumIsNotA);
	AllocIsNotA = NewIsNotSize;
	IdObject** tmp = new IdObjectPtr[AllocIsNotA];
	for (int i=NumIsNotA-1; i>=0; i--)
	    tmp[i] = IsNotAArray[i];
	delete [] IsNotAArray;
	IsNotAArray = tmp;
    }
}

bool IdObject::is_a(IdObject *x)
{
    if (x==this) return 1;
    if (NumIsA==0 && NumIsNotA==0 && InitIsA)
    {   // first call to is_a for this class
	(*InitIsA)(this);
    }
    // check to see if relationship known
    int i;
    for (i=NumIsA-1; i>=0; i--)
        if (x==IsAArray[i]) return 1;
    for (i=NumIsNotA-1; i>=0; i--)
        if (x==IsNotAArray[i]) return 0;
    
    // not known, do recursive search
    for (i=NumIsA-1; i>=0; i--)
	if (IsAArray[i]->is_a(x))
	{   add_is_a(x);
	    return 1;
	}
    
    // not found, record that this is not an x
    if (NumIsNotA >= AllocIsNotA) Realloc(AllocIsA, 2*AllocIsNotA+15);
    IsNotAArray[NumIsNotA++] =x;
    return 0;
}

void IdObject::add_is_a(IdObject *x)
{
    // check for and avoid redundant addition
    if (x==this) return;
    int i;
    for (i=NumIsA-1; i>=0; i--)
        if (x==IsAArray[i]) return;
    
    // check for inconsistency
    for (i=NumIsNotA-1; i>=0; i--)
        if (x==IsNotAArray[i])
	{   cerr << "Error: " << name() 
		<< " already queried and determined not to be an " 
		<< x->name()
		<< " and so can't do add_is_a\n" << flush;
	    abort();
	}

    if (NumIsA >= AllocIsA) Realloc(2*AllocIsA+15, AllocIsNotA);
    IsAArray[NumIsA++] =x;
}

NamedClass* NamedClass::read_new(istream &in)
{
    char word[1000];
    get_word(in, word,'=');
    if (!in.good() && word[0]==0)
	return 0;	// end of file or other lack of data, not an error
    
    if (!EqualStrings(word, "ClassName",1))
    {    cerr << "Error: need \"ClassName =\" followed by "
		<< "a registered class name\n"
		<<" to start reading a new object, but saw " 
		<< word << " instead.\n"
		<< flush;
	return 0;
    }

    return read_new_after_classname(in);
}

NamedClass * NamedClass::read_new_after_classname(istream &in)
{
    char word[1000];
    get_word(in, word);
    const IdObject *id = IdObject::id(word);
    if (!id)
    {    cerr << "Error: " << word << " is not a registered type.\n"
		<< "Can't create new objects of unregistered types.\n"
		<< flush;
	return 0;
    }
    NamedClass *read_into = id -> create();

    if (read_into && read_into->read_knowing_type(in))
	return read_into;
    
    get_word(in, word, '=');
    if (!EqualStrings(word, "EndClassName",1))
    {    cerr << "Warning: need \"EndClassName = " << id->name()
		<< "\" to end reading, but saw " << word << " instead.\n"
		<< flush;
	return read_into;
    }
    
    get_word(in,word);
    const IdObject *end_id = IdObject::id(word);
    if (end_id !=id)
    {    cerr << "Warning: EndClassName " << word << " doesn't match "
		<< id->name() << "\n" << flush;
    }
    return read_into;

}

// read many objects from istream (until EOF)
//	and pushes newly allocated objects on the end of the vector
// returns number read.
// If OnlyThisType is set, an objects that is not of the
//	specified type will cause and error.
int NamedClass::read_many_new(
	vector<NamedClass*>& objects,
	istream &named_in,
	IdObject* OnlyThisType)
{
    
    NamedClass *nc = NULL;
    int num_read = 0;

    for(nc = NamedClass::read_new(named_in);
    	named_in.good();
	nc = NamedClass::read_new(named_in))
    {   if (nc != NULL)
        {   if (OnlyThisType && !nc->is_a(OnlyThisType))
            {	cerr << "ERROR: NamedClass::read_many_new():"
                     << " required to read " << OnlyThisType->name()
		     << " , but got "
                     << nc->type()->name()
                     << "\n";
                abort();
            }
	    objects.push_back(nc);
            num_read++;
        }
    }
    return num_read;
}


// opens filename for input
//	and does read_man_new on the input stream
int NamedClass::read_many_new(vector<NamedClass*>& objects,
	const char *filename,
	IdObject* OnlyThisType)
{
    //cerr << "Entering NamedClass::read_many_new\n";

    gzifstream *named_in = Filenames::open_input(filename);
    if (!named_in) 
    {   cerr << "Error: can't open file for read_many_new(\"" << filename 
    		<< "\", " << OnlyThisType->name()
		<< "\n";
        return 0;
    }
       
    int num_read=read_many_new(objects, *named_in, OnlyThisType);
    named_in->close();
    delete named_in;

    return num_read;
}

istream& NamedClass::read(istream &in) 
{
    char word[1000];
    get_word(in, word,'=');
    if (!in.good() && word[0]==0)
	return in;	// end of file or other lack of data, not an error
    
    if (!EqualStrings(word, "ClassName",1))
    {    cerr << "Error: need \"ClassName =\" followed by "
		<<  type()->name()
		<<" to read an already allocated object of this type."
		<< "Input has '"
		<< word << "' instead of 'ClassName'.\n"
		<< flush;
	return in;
    }

    get_word(in, word);
    const IdObject *id = IdObject::id(word);
    if (!id || id!=type())
    {    cerr << "Error: trying to read class " << word 
		<< " into object of type " << type()->name()
		<< "\n"	<< flush;
	return in;
    }


    if (read_knowing_type(in)) return in;
    
    get_word(in, word, '=');
    if (!EqualStrings(word, "EndClassName",1))
    {    cerr << "Warning: need \"EndClassName = " << type()->name()
		<< "\" to end reading, but saw " << word << " instead.\n"
		<< flush;
    }
    get_word(in, word);
    id = IdObject::id(word);
    if (!id || id!=type())
    {    cerr << "Error: EndClassName= " << word 
		<< " Does not match ClassName= " << type()->name()
		<< "\n"	<< flush;
	return in;
    }

    
    return in;
}

const int BUF_LEN = 1023;


// skip to end of line
int ReadComment(istream &in,
		    NamedClass *change,
		    NamedClassCommand *self)
{   SkipSeparators(in, 1, '\n');
    return 1;
}

int ReadName(istream &in,
		    NamedClass *change,
		    NamedClassCommand *self)
{
    char word[BUF_LEN];
    get_word(in, word);
    NamedObject *nm_change = dynamic_cast<NamedObject*>(change);
    assert(nm_change);
    nm_change->set_name(word);
    return 1;
}

int VerifyClassName(istream &in, 
	    NamedClass *change,
	    NamedClassCommand* self)
{
    char word[BUF_LEN];
    get_word(in, word);
    const IdObject *end_id = IdObject::id(word);
    if (end_id != change->type())
    {    cerr << "Warning: " << self->name() << word << " doesn't match "
                << change->type()->name() << "\n" << flush;
    }
    // continue if "ClassName", stop if "EndClassName"
    return EqualStrings(self->name(), "ClassName", 1);
}


//CHANGE LOG:
// 3 Dec 1995 Kevin Karplus
//	added test for InitIsA before calling it in is_a()
// 4 Dec 1995 Kevin Karplus
//	added help string to constructor
// 5 Nov 1999 Kevin Karplus
//	split read_new_after_classname off from read_new
// Sat Jun 18 10:20:49 PDT 2005 Kevin Karplus
//	added read_many_new() functions
// Thu Jul  7 04:15:56 PDT 2005 Kevin Karplus
//	Missing OnlyThisType condition added to error-check in
//	read_many_new() 
// Sun Nov 26 09:21:30 PST 2006 Kevin Karplus
//	Changed type() to type()->name() in two error messages
// Sun Apr 19 16:48:54 PDT 2009 Kevin Karplus
//	Added NamedClassCommand functions ReadComment, ReadName,VerifyClassName
