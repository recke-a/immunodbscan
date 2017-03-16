/*
 * IdenticalSequencesMask.hpp
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


#include "BitsetVector.hpp"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/replace.h>
#include <thrust/copy.h>
#include <vector>
#include <string>
#include "StridedStringsOperations.hpp"


class device_data_container
{
	thrust::device_vector<unsigned int> reference_vector;
	thrust::device_vector<unsigned int> n_copies_vector;
	
	bitset_vector     	Vgenes;	
	bitset_vector 		Dgenes;
	bitset_vector 		Jgenes;
	
	StridedStringSet 	SequencesSet;
		
	public:
	device_data_container(std::size_t length, bitset_vector&& _Vgenes, bitset_vector&& D_genes, bitset_vector&& _Jgenes, std::vector<std::string>&& inputStrings);
	
	void mask_by_query_index(unsigned int query_index);
		
	void get_3_epsilon_approximate_sequences_with_mask( thrust::device_vector<unsigned int>::iterator currentLeader, 
	                                                    thrust::device_vector<unsigned int>& LeaderIndices, thrust::device_vector<bool>&  MaskVector); 
	                                                    
	template<typename bitset_iterator_type, typename widths_type> void get_2_epsilon_approximate_sequences_with_mask(
		thrust::device_vector<unsigned int>::iterator currentLeader, 
		thrust::device_vector<unsigned int>& LeaderIndices, 
		thrust::device_vector<bool>& MaskVector, 
		thrust::device_vector<unsigned int>& leader_1_indices, 
		unsigned int N_indices, 
		bitset_iterator_type GeneMasksContainer_iterator,
		widths_type GeneMasksContainer_widths);
		
};

