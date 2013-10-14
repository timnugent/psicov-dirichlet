// OrderSort.cc
// Kevin Karplus
// 9 June 1995

// ABSTRACT
// OrderSort takes an array of floats (value) and
//	returns a newly allocated array of integers (index)
//	such that value[index[i+1]] >= value[index[i]]

#include "Prob.h"

// Sort an index array (containing the values 0..num_to_sort-1)
// so that value[index[i]] is in ascending order.
// This sort is a simple insertion sort.
// Sort is stable (stability guranteed by extra test when keys are equal).
template <typename T>
void InsertSortIndex(int *index, const T * value, int num_to_sort)
{
    for (int i=1; i<num_to_sort; i++)
    {   int move=index[i];
        T inserting=value[move];
	int j;
	for(j=i-1; 
		j>=0 && (value[index[j]]>inserting ||
			 (value[index[j]]==inserting && index[j]>move));
		j--)
	{    index[j+1] = index[j];
	}
	index[j+1] = move;
    }
}

const int MAX_LENGTH_FOR_INSERT_SORT = 40;

template <typename T>
void SortIndex(int *index, const T * value, int num_to_sort)
{
    int move;	// temporary for swapping index values
    switch(num_to_sort)
    {    case 0: return;
         case 1: return;
         case 2: if (value[index[0]] > value[index[1]]
		 	|| (value[index[0]] == value[index[1]]
				&& index[0] >index[1]))
	 	{    move=index[0];
		    index[0]=index[1];
		    index[1]=move;
		}
		return;
    	case 3:  if (value[index[0]] > value[index[1]]
		 	|| (value[index[0]] == value[index[1]]
				&& index[0] >index[1]))
	 	{    move=index[0];
		    index[0]=index[1];
		    index[1]=move;
		}
		if (value[index[1]] > value[index[2]]
			|| (value[index[1]] == value[index[2]]
				&& index[1]>index[2] ))
	 	{    move=index[2];
		    index[2]=index[1];
		    index[1]=move;
		    if (value[index[0]] > value[index[1]]
		 	|| (value[index[0]] == value[index[1]]
				&& index[0] >index[1]))
	 	    {   move=index[0];
		    	index[0]=index[1];
		    	index[1]=move;
		    }
		}
		return;	
    }
    if (num_to_sort<MAX_LENGTH_FOR_INSERT_SORT) 
    {   InsertSortIndex<T>(index,value,num_to_sort);
    	return;
    }

    int pivot_loc =num_to_sort/2;
    T pivot= value[index[pivot_loc]];
    
    move = index[pivot_loc];
    index[pivot_loc] = index[num_to_sort-1];	
    index[num_to_sort-1] = move;	// guarantee a stopping point
    
    int left=0; 
    int right=num_to_sort-2;
    while(left<=right)
    {   while( value[index[left]]<pivot 
    		|| (value[index[left]]==pivot && index[left]<move))
	    left++;
//        assert (left <num_to_sort);
	while ( (value[index[right]] > pivot 
		    || (value[index[right]]==pivot && index[right]>=move))
		&& right>left)
	    right--;
	if (left>=right) break;
	int tmp=index[left];
	index[left]=index[right];
	index[right]=tmp;
	left++;
	right--;
	if (left>right) break;
    }
    
    move = index[num_to_sort-1];
    index[num_to_sort-1] = index[left];
    index[left] = move;
//    assert(value[index[num_to_sort-1]] >= pivot);
    
    int num_before_pivot = left;
    SortIndex(index,value,num_before_pivot);
    SortIndex(index+num_before_pivot+1,value, num_to_sort-num_before_pivot-1);
}

template <typename T>
int * CreateAndSortIndex(const T * value, int num_to_sort)
{
    if (num_to_sort <= 0)
    	return NULL;

    // create the index array
    int * index= new int[num_to_sort];
    for (int i=num_to_sort-1; i>=0; i--)
    	index[i] = i;

    SortIndex<T>(index,value,num_to_sort);
    return index;
}

int * OrderSort(const int * value, int num_to_sort)
{
    return CreateAndSortIndex<int>(value,num_to_sort);
}

int * OrderSort(const float * value, int num_to_sort)
{
    return CreateAndSortIndex<float>(value,num_to_sort);
}

int * OrderSort(const double * value, int num_to_sort)
{
    return CreateAndSortIndex<double>(value,num_to_sort);
}

int * OrderSort(const Prob * value, int num_to_sort)
{
    return CreateAndSortIndex<Prob>(value,num_to_sort);
}

// CHANGE LOG:
// 10 April 2004 Kevin Karplus
//	OrderSort rewritten to use template-version of quicksort,
//	rather than calling qsort.
//	Note: template usage hidden in .cc file, not visible in .h file.
// 27 April 2004 Kevin Karplus
//	modified tests to make OrderSort stable.
// 28 May 2004 Kevin Karplus
//	Added declared, but missing, OrderSort(const int*, int)
