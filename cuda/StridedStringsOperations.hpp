/*
 * StridedStringsOperations.hpp
 * 
 * Copyright 2017 Andreas Recke <andreas@AntigoneLinux>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

// This program part defines a strided stringset and the kernel for comparison as well as translation of DNA into amino acid codes

#pragma once
#include <thrust/device_vector.h>
#include <vector>
#include <thrust/host_vector>
#include <string>





class StridedStringSet
{
	// this class stores and works on 
	
	thrust::device_vector<std::size_t> StringHash; // may have an interesting use
	thrust::device_vector<char>  Strings;    // the data 
	thrust::device_vector<unsigned int> StringLengths;
	float threshold;
	bool use_Hamming;

	unsigned n_entries;
	unsigned int width;            // the 
	
	StridedStringSet(std::vector<std::string>&& inputStrings, float _threshold);
	
	void get_equal_strings(thrust::device_vector<bool>& result, unsigned int query_index);
	
	void get_3_epsilon_approximate_sequences_with_mask( thrust::device_vector<unsigned int>::iterator currentLeader, 
                                                    thrust::device_vector<unsigned int>& LeaderIndices, 
                                                    thrust::device_vector<bool>& MaskVector);
	
};
