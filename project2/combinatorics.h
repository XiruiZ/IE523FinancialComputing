/*
 *  combinatorics.h
 *  Created by Prof. Ramavarapu Sreenivas 
 *
 */

#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace std;

template <typename T>
vector < vector < T > > Pick_R(vector < T > array, int R)
{
	// this is the final result
	vector <vector <T> > result;
	
	if ((R <= array.size()) && (R > 0))
	{
		// Create a boolean array where the first r items are "true" (i.e. 1)
		// the remaining n-r items are "false" (i.e. 1)
		vector <bool> did_i_pick_this(array.size());
		fill(did_i_pick_this.begin(), did_i_pick_this.begin() + R, true);
		
		// I am going to use the prev_permutation STL in C++ (see 
		// http://www.cplusplus.com/reference/algorithm/prev_permutation/ )
		// to generate the required combinatorial outputs 
		
		do 
		{
			vector <T> current_set;
			for (int i = 0; i < array.size(); ++i) 
			{
				if (did_i_pick_this[i]) 
					current_set.push_back(array[i]);
			}
			result.push_back(current_set);
		} while (prev_permutation(did_i_pick_this.begin(), did_i_pick_this.end())); 
	}
	
	return (result);
}

// print the vector-of-vectors that is output by the selection process
template <typename TYPE>
void print_a_vector_of_vectors(vector < vector <TYPE> > array)
{
	
	for (int i = 0; i < array.size(); i++) {
		cout << setw(3) << i+1 << ": { ";
		for (int j = 0; j < array[i].size(); j++) 
			cout << array[i][j] << " ";
		cout << "}" << endl;
	}
}

// Power set -- computes the non-empty subsets of an arbitrary set
template <typename T>
vector < vector < T > > Power_Set(vector < T > array)
{
	// this is the final result
	vector <vector <T> > result;
	
	for (int R = 1; R <= array.size(); R++)	
	{
		// temporary
		vector <vector <T> > result1;
		result1 = Pick_R(array, R);
		for (int i = 0; i < result1.size(); i++)
			result.push_back(result1[i]);
	}
	return (result);
}

#endif