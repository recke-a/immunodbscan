/*
 * LeaderImplementation.hpp
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
#include "BitsetVector.hpp"

class FMDBSCANclass
{
	// This is leader data
	thrust::device_vector<unsigned int> LeaderIndices;            // a list of indices pointing into the sequence list
	thrust::device_vector<bool> UnVisitedLeaderTags;        // corresponding to the List of leaders for clustering, will be true if a leader hasn't been clustered yet
	
	bitset_vector    VgeneLeaderMasks;  // this is for candidate searches! Not to prove direct density reachability between two leader regions
	bitset_vector    JgeneLeaderMasks;  // this is for candidate searches!
		
	// The following data is the original data
	device_data_container sequencing_data;
	thrust::device_vector<int> LeaderAssociations;    // a vector containing leader information in the total data set
	thrust::device_vector<int> LeaderSizes;    // a vector containing leader information in the total data set
	
	// This is the annotation data
	thrust::device_vector<int>          is_core_point_property; 
	thrust::device_vector<int>          ClusterIds;                // the final list of ClusterIds
		

	    
	// This is temporary data ... to reduce new and free operations
	thrust::device_vector<unsigned int> leader_1_indices;
	thrust::device_vector<unsigned int> leader_2_indices;
	
	thrust::device_vector<unsigned int> candidates;
	thrust::device_vector<bool> candidates_mask;
	
	
	// further data
	bool use_J_gene;
	float threshold;
	unsigned int currentLeader;
	unsigned int N_sequences;
	unsigned int N_leaders;
	
	public:
	
	
	
	
	
	
	
	
	struct equal_to_value_op : public thrust::unary_function<unsigned int, bool>
	{
		private:
		unsigned int comp_value;
		
		public:
		equal_to_value_op(unsigned int _comp) : comp_value(_comp) {}
		
		bool operator()(unsigned int arg)
		{
			return arg == comp_value;
		}
		
	};
	
	typedef thrust::tuple<bool, unsigned int> unvisited_core_t;
	
	struct is_unvisited_core_op : public thrust::unary_function<unvisited_core_t, bool>
	{
		bool operator()(unvisited_core_t arg)
		{
			return arg.get<0>() && ( arg.get<1> == 1u );
		}
	};
	
	
	typedef thrust::tuple<bool, bool> two_bool_values_t;
	struct zipped_logical_and : public thrust::unary_function<two_bool_values_t, bool>
	{
		bool operator()(two_bool_values_t arg)
		{
			return arg.get<0>() && arg.get<1>();
		}
	}; 
	
	
};

//*********************************** Implementations *********************************




/*
 * 
 * name: get_3_epsilon_environment
 * @param: none
 * @return: none
 * 
 * This function selects candidate leaders that may be connected with the "currentLeader" from the list of all unvisited leaders.
 * It modifies the class member "MaskVector" that is a logical vector which is either true or false for a certain leader
 * "currentLeader" is the index of the leader pointing into the list of "LeaderIndices" that then points to the original sequencing data
 */
void FMDBSCANclass::get_3_epsilon_environment()
{
	using namespace thrust;
		
	// the first of the candidates is always the leader from which we are querying!!!	
		
	// 1: get the leaders that are still not visited
	auto ende = copy_if(device, make_counting_iterator(0u), make_counting_iterator(0u) + N_leaders, UnVisitedLeaderTags.begin(), candidates.begin(), identity<bool>());
	
	// 2: determine number of leaders in this list
	unsigned int n_candidates = ende - candidates.begin();
	
	// 3: generate an iterator that points to the strings were interested in
	auto leader_index_to_seq_index_iterator = make_permutation_iterator(LeaderIndices.begin(), candidates.begin());
	auto str_index_iterator = make_permutation_iterator(sequencing_data.string_iterator_begin(), leader_index_to_seq_index_iterator);
	
	// sequencing_data.string_iterator_begin() is a jumping iterator combined with a string lengths iterator
	
	epsilon_3_functor distance_functor(threshold);
	
	// 4: prepare mask vector
	if (use_J_gene)
	{
		auto gene_index_iterator = make_permutation_iterator(make_zip_iterator(make_tuple(VgeneLeaderMasks.begin(), JgeneLeaderMasks.begin())), candidates.begin());
				
		call_query_kernel(gene_index_iterator, str_index_iterator, n_candidates, candidates_mask.begin(), distance_functor);
	} 
	else {
		auto gene_index_iterator = make_permutation_iterator(make_zip_iterator(make_tuple(VgeneLeaderMasks.begin())), candidates.begin());
				
		call_query_kernel(gene_index_iterator, str_index_iterator, n_candidates, candidates_mask.begin(), distance_functor);
	}
	
	candidates.erase(device, candidates.remove_if(device, candidates.begin(), candidates.end(), candidates_mask.begin(), identity<bool>()), candidates.end());
	
	// ready!
}
	
/*	
	oder soo:::
	
	
	// 1st step: the list of unvisited leaders are copied into the mask
	thrust::copy(thrust::device, UnVisitedLeaderTags.begin(), UnVisitedLeaderTags.end(), MaskVector.begin());
	
	// 2nd step: the mask is reduced to contain only leaders that match the V (and J) gene distribution
	// Important here is that the LeaderMasks contain all possible combinations of V and J gene bits of the leader region members
	VgeneLeaderMasks.get_equivalent_entries_with_mask(currentLeader, MaskVector);
	if (use_J_gene) JgeneLeaderMasks.get_equivalent_entries_with_mask(currentLeader, MaskVector);
			
	// 3rd step: select the Leaders that may be connected with two intermediate elements (3 epsilon environment)
	sequencing_data.get_3_epsilon_approximate_sequences_with_mask( (LeaderIndices.begin() + currentLeader), LeaderIndices, MaskVector);  // searches within Leaders!!! 
	
	// Completed
}

*/




void FMDBSCANclass::get_2_epsilon_environment()
{
	// we're searching for indices pointing to the sequence data where the leader number fits and the sequence is a core point
	equal_to_value_op comparison_to(currentLeader);
	
	// the following fance iterator helps to fuse GPU kernels
	auto interesting_leader_1_indices_iter = thrust::make_transform_iterator(
												thrust::make_zip_iterator(
													thrust::make_tuple(
														thrust::make_transform_iterator(LeaderSizesAndAssociations.begin(), comparison_to),
														is_core_point_property.begin()
													)
												),
												zipped_logical_and()
											 );
		
	// this is a simple trick to translate the logical vector to a list of indices
	auto ende_target_indices = thrust::copy_if(thrust::device,
											   thrust::make_counting_iterator(0), 
											   thrust::make_counting_iterator(0) + N_sequences,
											   interesting_leader_1_indices_iter,
											   leader_1_indices.begin(),
											   thrust::identity<bool>());
											   
    unsigned int N_indices = ende_target_indices - leader_1_indices.begin();
    
    // the following order is to mask all leader regions that cannot be connected to the current leader reagion
    if (use_J_gene) 
    {
		sequencing_data.get_2_epsilon_approximate_sequences_with_mask(currentLeader, LeaderIndices, MaskVector, leader_1_indices, N_indices, 
			thrust::make_zip_iterator(thrust::make_tuple(VgeneLeaderMasks.container.begin(), JgeneLeaderMasks.container.begin())), 
			thrust::make_tuple(VgeneLeaderMasks.width, JgeneLeaderMasks.width));		
	} else
	{
		sequencing_data.get_2_epsilon_approximate_sequences_with_mask(currentLeader, LeaderIndices, MaskVector, leader_1_indices, N_indices, 
			VgeneLeaderMasks.container.begin(),
			VgeneLeaderMasks.width);
	}
}







void FMDBSCANclass::expand_cluster(unsigned int query_index)
{
	// set current leader
	currentLeader = query_index;
	
	// unmark this leader as unvisited
	UnVisitedLeaderTags[currentLeader] = false;
	
	// get candidates
	get_3_epsilon_environment();

	// further reduce candidates
	get_2_epsilon_environment();
	
	
	
	
	
}


void FMDBSCANclass::next_cluster()
{
	// each leader must be a core point
	// the core point property has been set for all leaders as 1 or 0. Otherwise this property is -1 (for undefined)
	// search for next leader that is not part of any cluster and that is a core point
		
	auto core_perm_it = thrust::make_permutation_iterator(is_core_point_property.begin(), LeaderIndices.begin());
	auto visit_core_zip_iter = thrust::make_zip_iterator(thrust::make_tuple(UnVisitedLeaderTags.begin(), core_perm_it);
	
	unsigned int query_index = thrust::find_if(thrust::device, visit_core_zip_iter, visit_core_zip_iter + N_sequences, is_unvisited_core_op()) - visit_core_zip_iter;
	
	expand_cluster(query_index);
}
