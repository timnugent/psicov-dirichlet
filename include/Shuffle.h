#ifndef _SHUFFLE_H_
#define _SHUFFLE_H_

#include "Random.h"

// two tiny routines for creating a random permutation.

// randomly shuffle an array of num ints
inline void shuffle(int *cards, int num)
{
    for (int r=num-1; r>=0; r--)
    {   int s = irandom(num);
    	assert (0<=s && s<num);
	int tmp = cards[r];
	cards[r] = cards[s];
	cards[s] = tmp;
    }
}

// create an array with a random permutation of 0..num-1
inline int * random_permute(int num)
{
    assert(num>0);
    int* order= new int[num];	// randomized order for rotamers
    for (int r=num-1; r>=0; r--)
       	order[r] = r;
    shuffle(order, num);
    return order;
}



// CHANGE LOG:
// Sat Jan 15 12:28:58 PST 2005 Kevin Karplus
// 	added include for Random.h

#endif
