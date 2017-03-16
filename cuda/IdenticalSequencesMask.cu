/*
 * IdenticalSequencesMask.cu
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

#include "IdenticalSequencesMask.hpp"

/************ Implementations ****************************/

device_data_container::device_data_container(std::size_t length, bitset_vector&& _Vgenes, bitset_vector&& D_genes, bitset_vector&& _Jgenes, std::vector<std::string>&& inputStrings) :
	Vgenes(_Vgenes),
	Dgenes(_Dgenes),
	Jgenes(_Jgenes),
	SequencesSet(std::forward(inputStrings)),
	n_copies_vector(_Vgenes.n_elem, 1u),
	reference_vector(_Vgenes.n_elem)   // init reference_vector with 1 elements
{
	auto counter = thrust::make_counting_iterator<unsigned int>(0);
	thrust::copy(counter, counter + _Vgenes.n_elem, reference_vector.begin());
	
	auto n_copies_vector_iterator = n_copies_vector.begin();
	
	while (n_copies_vector_iterator != n_copies_vector.end())
	{
		n_copies_vector_iterator = thrust::find(thrust::device, n_copies_vector_iterator, n_copies_vector.end(), 1u);  // search for next position of vector
		if (n_copies_vector_iterator != n_copies_vector.end())
		{
			unsigned int query_index = n_copies_vector_iterator - n_copies_vector.begin();
			mask_by_query_index(query_index);
		}
	}
}


void device_data_container::mask_by_query_index(unsigned int query_index)
{
	thrust::device_vector<bool> equivalents(reference_vector.size(), true);
	Vgenes.get_identical_entries_with_mask(query_index, equivalents);
	Dgenes.get_identical_entries_with_mask(query_index, equivalents);
	Jgenes.get_identical_entries_with_mask(query_index, equivalents);
	SequencesSet.get_equal_strings(query_index, equivalents);
	
	unsigned int n_counts = thrust::count_if(thrust::device, equivalents.begin(), equivalents.end(), thrust::identity<bool>());
	thrust::replace_if(thrust::device, reference_vector.begin(), reference_vector.end(), equivalents.begin(), thrust::identity<bool>(), query_index);
	thrust::replace_if(thrust::device, n_copies_vector.begin(), n_copies_vector.end(), equivalents.begin(), thrust::identity<bool>(), 0u);
	n_copies_vector[query_index] = n_counts;
}


void device_data_container::get_3_epsilon_approximate_sequences_with_mask( thrust::device_vector<unsigned int>::iterator currentLeader, 
                                                    thrust::device_vector<unsigned int>& LeaderIndices, 
                                                    thrust::device_vector<bool>&  MaskVector)
{
	SequencesSet(currentLeader, LeaderIndices, MaskVector);
}



template<typename bitset_iterator_type, typename widths_type> void device_data_container::get_2_epsilon_approximate_sequences_with_mask(
		thrust::device_vector<unsigned int>::iterator currentLeader, 
		thrust::device_vector<unsigned int>& LeaderIndices, 
		thrust::device_vector<bool>& MaskVector, 
		thrust::device_vector<unsigned int>& leader_1_indices, 
		unsigned int N_indices, 
		bitset_iterator_type GeneMasksContainer_iterator,
		widths_type GeneMasksContainer_widths)
{
	
}
