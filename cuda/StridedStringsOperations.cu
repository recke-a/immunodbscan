/*
 * StridedStringsOperations.cu
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

#include "StridedStringsOperations.hpp"
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>


struct str_len_functor : public thrust::unary_function< std::string, unsigned int>
{
	__host__ unsigned int operator()(std::string x)
	{
		return x.size();
	}
};



// take up data

StridedStringSet::StridedStringSet(std::vector<std::string>&& inputStrings) :
  StringHash(inputStrings.size()),
  StringLengths(inputStrings.size())
{
	// translate Strings to hashes
	thrust::host_vector<std::size_t> temp_hash(inputStrings.size());
	thrust::transform(inputStrings.begin(), inputStrings.end(), temp_hash.begin(), std::hash<std::string>());
	thrust::copy(temp_hash.begin(), temp_hash.begin()+inputStrings.size(), StringHash.begin()); 
	
	thrust::host_vector<unsigned int> temp_lens(inputStrings.size());
	thrust::transform(inputStrings.begin(), inputStrings.end(), temp_lens.begin(), str_len_functor());
	thrust::copy(temp_lens.begin(), temp_lens.begin()+inputStrings.size(), StringLengths.begin()); 
	
	// find maximum string length
	width = thrust::reduce(thrust::device, StringLengths.begin(), StringLengths.end(), thrust::maximum<unsigned int>());
	
	// allocate Container and init with zero
	thrust::host_vector<char> temp_strings(inputStrings.size() * width, 0);
	
	
	// and copy data
	unsigned position = 0;
	for (auto s_it = inputStrings.begin(); s_it != inputStrings.end(); ++s_it, position += width) 
	    thrust::copy(s_it->begin(), s_it->end(), temp_strings.begin() + position);
    
    // Strings    = thrust::device_vector<char> (inputStrings.size() * width);
    Strings    = thrust::device_vector<char> (temp_strings.begin(), temp_strings.end());  // mal sehen!
    // thrust::copy(temp_strings.begin(), temp_strings.end(), Strings.begin());	
}




//******************************************************  A function that filters equal strings ****************************************************************************************

__global__ void equality_kernel(
								unsigned int query_index, 
								thrust::device_vector<char>::iterator str_subj_base, 
								thrust::device_vector<unsigned int>::iterator str_len_base,
								unsigned int string_width,
								thrust::device_vector<std::size_t>::iterator hash_base,
								thrust::device_vector<bool>::iterator result_base
							)
{
	unsigned int entry         = blockIdx.x;
		
	bool* result_pointer = thrust::raw_pointer_cast( & (*(result + entry)));
	
	if (*result_pointer) 
	{
		if ( *(hash_base + *(index_base + entry)) == *(str_subj_base + query_index) )
		{ 
			unsigned int strlengthQuery   = *(str_len_base + query_index);
			unsigned int strlengthSubject = *(str_len_base + *(index_base + entry));
			
			if (strlengthQuery != strlengthSubject)
			{
				*result_pointer = false;
			}
			else
			{
				int hamming;
				
				bool neq = false;
				if (threadIdx.x < strlengthQuery) 
				    neq = *(str_subj_base + ((*(index_base + entry))*string_width) + threadIdx.x) != *(str_subj_base + ((*(index_base + query_index))*string_width)  + threadIdx.x); 
				
				hamming = __synchthreads_count(neq);
				
				if (threadIdx.x == 0) *result_pointer = hamming == 0;
			}
		}	
	}
}



void StridedStringSet::get_equal_strings(thrust::device_vector<bool>& result, unsigned int query_index)
{
	equality_kernel<<<Strings.size(), width>>>(query_index, Strings.begin(), StringLengths.begin(), width, StringHash.begin(), result.begin());
}



//******************************************************  A function that filters the 3 eps environment ****************************************************************************************

template <typename T> __host__ __device__ T maximum(T a, T b)
{
	return (a>b)? a : b;
}

template <typename T> __host__ __device__ T minimum(T a, T b)
{
	return (a<b)? a : b;
}



struct epsilon_1_functor
{
	float threshold;
	
	epsilon_1_functor(float _threshold) : threshold(_threshold) {}
	
	__host__ __device__ unsigned int operator()(unsigned int query_len, unsigned int subj_len)
	{
		float fschwelle = maximum(float(query_len), float(subj_len)) * threshold;
		return fschwelle;
	}
};


struct epsilon_2_functor
{
	float threshold;
	
	epsilon_2_functor(float _threshold) : threshold(_threshold) {}
	
	__host__ __device__ unsigned int operator()(unsigned int query_len, unsigned int subj_len)
	{
		float inter_b_len = minimum(float(query_len)/(1.0f - threshold), float(subj_len)/(1.0f - threshold));
		float fschwelle = ( maximum(float(query_len), inter_b_len) + maximum(inter_b_len, float(subj_len)) ) * threshold;
		return fschwelle;
	}
};


struct epsilon_3_functor
{
	float threshold;
	
	epsilon_3_functor(float _threshold) : threshold(_threshold) {}
	
	__host__ __device__ unsigned int operator()(unsigned int query_len, unsigned int subj_len)
	{
		float inter_b_len = minimum(float(query_len)/(1.0f - threshold), float(subj_len)/(1.0f - threshold)/(1.0f - threshold));
        float inter_c_len = minimum(float(subj_len)/(1.0f - threshold), float(query_len)/(1.0f - threshold)/(1.0f - threshold));
			            			            
		float fschwelle = (maximum(float(query_len), inter_b_len) + maximum(inter_c_len, inter_b_len) + maximum(float(subj_len), inter_c_len))*threshold;
		return fschwelle;
	}
};


//******************************************************************************************************************

// first a simple hamming kernel
template <typename threshold_functor_type> __global__ void hamming_kernel(  
										thrust::device_vector<unsigned int>::iterator query_address,
										thrust::device_vector<unsigned int>::iterator index_base,
										thrust::device_vector<char>::iterator str_subj_base, 
										thrust::device_vector<unsigned int>::iterator str_len_base,
										unsigned int string_width,
										thrust::device_vector<bool>::iterator result,   // the result is initialized with true
										threshold_functor_type epsilon_functor   
                                     )
{
	unsigned int entry         = blockIdx.x;
	
	// copy pointer to the result from global into local memory for faster access
	bool* result_pointer = thrust::raw_pointer_cast( & (*(result + entry)));
	unsigned int query_index = *query_address;
	unsigned int subj_index  = *(index_base + entry);
	
				
	if (*result_pointer)
	{ 
		if ( *(hash_base + subj_index) != *(str_subj_base + query_index) )
		{
			unsigned int strlengthQuery   = *(str_len_base + query_index);
			unsigned int strlengthSubject = *(str_len_base + subj_index);
				
			if (strlengthQuery != strlengthSubject)
			{
				*result_pointer = false;
			} 
			else
			{
				int hamming;
				int schwelle  = epsilon_functor(strlengthQuery, strlengthSubject);  // must always be called with two values
					
				bool neq = false;
				if (threadIdx.x < strlengthQuery) 
				    neq = *(str_subj_base + (subj_index*string_width + threadIdx.x)) != *(str_subj_base + (query_index*string_width + threadIdx.x)); 
						
				hamming = __synchthreads_count(neq);
					
				if (threadIdx.x == 0) *result_pointer = hamming <= schwelle; 
			}	
		}
	}
}



// some important issues:
// 1. the leaders should all be different from each other, so Hash values are unimportant, but the Hamming distance should work
// 2. the function to calculate the threshold is a bit difficult
// 3. we use a one-dimensional grid and block setting

// the following kernel must be tested separately!!

template <typename threshold_functor_type, bool onlyHamming=false> __global__ void index_based_within_epsilon_kernel (  
										thrust::device_vector<unsigned int>::iterator query_address,
										thrust::device_vector<unsigned int>::iterator index_base,
										thrust::device_vector<char>::iterator str_subj_base, 
										thrust::device_vector<unsigned int>::iterator str_len_base,
										unsigned int string_width,
										thrust::device_vector<bool>::iterator result,
										threshold_functor_type epsilon_functor   
                                     )
{
	
	// request shared memory
	extern __shared__ char s[];
	
	unsigned int entry     = blockIdx.x;
	unsigned int position = threadIdx.x;
	bool* result_pointer = thrust::raw_pointer_cast( & (*(result + entry)));
	
	unsigned int query_index = *query_address;
	unsigned int subj_index  = *(index_base + entry);
	
	if (*result_pointer)
	{
		
		// distribute shared memory
		char* subjstring  = s;
		char* querystring = s[blockDim.x];
		unsigned short *line_values   = (unsigned short*)&querystring[blockDim.x];
		unsigned short *intermediate  = &line_values[blockDim.x];
		bool           *non_identity  = (bool*)&intermediate[blockDim.x];
		
		// local memory
		unsigned int position  = threadIdx.x;
		unsigned int query_len = *(str_len_base + query_index),
			          subj_len = *(str_len_base + subj_index);
			 
		unsigned int schwelle = epsilon_functor(query_len, subj_len);
		
		if (query_len + schwelle > subj_len || subj_len + schwelle > query_len)
		{
			if (position == 0) *result_pointer = false;
			return;
		}
		else
		{
			char* query_string_pointer = thrust::raw_pointer_cast( & (*(str_len_base + (query_index * string_width))));
			char* subj_string_pointer  = thrust::raw_pointer_cast( & (*(str_len_base + (subj_index  * string_width))));
			
			// copy strings into shared memory
			subjstring[position]  = subj_string_pointer[position];
			querystring[position] = query_string_pointer[position];
				
			// then calculate Hamming distance
										
			bool neq;
			if (position < query_len) neq = subjstring[position] != querystring[position]; else neq = false;
			unsigned int hamming = __syncthreads_count(neq);  // works only, if there is one block per string!!!
			
			if (onlyHamming)
			{
				*result_pointer = query_len == subj_len && hamming <= schwelle;
				return;
			} else 
			{
				if (query_len == subj_len && hamming <= schwelle) 
				{
					if (position == 0) *result_pointer = true;
				    return;
				}
				else
				{							
					// prepare 1st line
					line_values[position] = position;
					
					// walk throught the characters of the query string
					for (unsigned int str2pos = 0; str2pos < query_len; ++str2pos)   
					{
					    non_identity[position] = false;	
					    intermediate[position] = str2pos + 1;
					
						if (position > 0) 
						{
							if (subjstring[position-1] == querystring[str2pos]) {
								intermediate[position] = line_values[position-1];
							}
							else
							{
						        intermediate[position] = 1 + minimum(line_values[position-1], line_values[position]);
						        non_identity[position] = true;
							}					    
						}	    
						
						__syncthreads();  // finish the above stuff and go to the next step, logically applies to all threads!!!!
						    
						unsigned int i;    
						// this code uses a partitioning of the logical strucure in the Levenshtein distance matrix
						if (position <= subj_len) 
						{
						    if (!non_identity[position]) {
								i = position + 1;
								line_values[position] = intermediate[position];
											    
							    while (non_identity[i] && i <= subj_len)
							    {
									line_values[i] = minimum(intermediate[i], line_values[i-1] + 1);
									++i;
								}
							}
						}
						__syncthreads();  //  logically applies to all threads!!!!
							
					}
					// write result
					if (position == 0) *result_pointer = line_values[subj_len] <= schwelle;
				}	
			}
		}
	}
}

//**************************************************************************************************

void SequencesSet::get_3_epsilon_approximate_sequences_with_mask( thrust::device_vector<unsigned int>::iterator currentLeader, 
                                                    thrust::device_vector<unsigned int>& LeaderIndices, 
                                                    thrust::device_vector<bool>&  MaskVector)
{
	//
	using namespace thrust;
	unsigned int stringlength = 1 + width / 32;
	stringlength *= 32;
	
	unsigned int shared_memory_usage = (sizeof(char) + sizeof(char) + sizeof(unsigned short) + sizeof(unsigned short) + sizeof(bool)) * stringlength;
	
	epsilon_3_functor eps_fun(threshold);
	
	if (use_hamming) {
		//hamming_kernel<<<n_entries, width+1>>>(
			//currentLeader, LeaderIndices.begin(), 
		    //Strings.begin(), StringLengths.begin(), width, MaskVector.begin(), 
		    //eps_fun
		//);
		index_based_within_epsilon_kernel<epsilon_3_functor, true><<<n_entries, width+1, shared_memory_usage>>>(
			currentLeader, LeaderIndices.begin(), 
		    Strings.begin(), StringLengths.begin(), width, MaskVector.begin(), 
		    eps_fun
		);
	} else 
	{	
		index_based_within_epsilon_kernel<<<n_entries, width+1, shared_memory_usage>>>(
			currentLeader, LeaderIndices.begin(), 
		    Strings.begin(), StringLengths.begin(), width, MaskVector.begin(), 
		    eps_fun
		);
	}	
}


//**************************************************************************************************

template <unsigned int selector, typename iterator_tuple_type> __device__ 
bool impl_compare_bitset_values(iterator_tuple_type bitvalues1, iterator_tuple_type bitvalues2)
{
	if (selector == 1)
	{
		typedef thrust::iterator_traits<thrust::tuple_element<0, iterator_tuple_type> >::value_type iterator_value_type;   // for correct generic programming ... should be int!
		bool result = true;
		iterator_value_type *pointer1 = thrust::raw_pointer_cast( & (*(bitvalues1.get<0>())) );
		iterator_value_type *pointer2 = thrust::raw_pointer_cast( & (*(bitvalues2.get<0>())) );
		unsigned int width = bitvalues1.get<0>().get_jumping_length();
		
		for (unsigned int i = 0; i < width; ++i) result = result || ( pointer1[i] & pointer2[i])); 
		return result;
	} else {
		bool result = compare_bitset_values<selector-1>(bitvalues1, bitvalues2);
		if (result) 
		{
			typedef thrust::iterator_traits<thrust::tuple_element<selector-1, iterator_tuple_type> >::value_type iterator_value_type;
			iterator_value_type *pointer1 = thrust::raw_pointer_cast( & (*(bitvalues1.get<selector-1>())) );
			iterator_value_type *pointer2 = thrust::raw_pointer_cast( & (*(bitvalues2.get<selector-1>())) );
			unsigned int width = bitvalues1.get<selector-1>().get_jumping_length();
			
			for (unsigned int i = 0; i < width; ++i) result = result || ( pointer1[i] & pointer2[i])); 
		}
		return result;
	}
}

template <typename iterator_tuple_type> __device__ 
bool compare_bitset_values(iterator_tuple_type bitvalues1, iterator_tuple_type bitvalues2)
{
	return impl_compare_bitset_values< tuple_size<iterator_tuple_type>::value >(bitvalues1, bitvalues2);
}


// the following kernel uses a 2D grid with a 1D block to check whether the points in Leader region 1 may be
// connected with the points in Leader region 2.
// if a hit is found, this is checked.
// the result is a list of leader indices which may are checked 

template <bool onlyHamming=false, bool second_gene_entries_are_leaders=false, typename threshold_functor_type, typename my_iterator_type > 
__global__ 
void leader_density_reachable (  
										

										thrust::device_vector<unsigned int>::iterator index_1_base, 
										thrust::device_vector<unsigned int>::iterator index_2_base, 
										
										my_iterator_type genes_base_1,   // multiple jumping iterators!
										my_iterator_type genes_base_2,   // multiple jumping iterators!
										
										thrust::device_vector<char>::iterator str_subj_base, 
										thrust::device_vector<unsigned int>::iterator str_len_base,
										thrust::device_vector<unsigned int>::iterator leader_association_base,
										unsigned int string_width,
										thrust::device_vector<bool>::iterator result,
										threshold_functor_type epsilon_functor
                                     )
{
	
	// request shared memory
	extern __shared__ char s[];
	
	unsigned int entry1      = *(index_1_base + blockIdx.x);  // the query
	unsigned int entry2      = *(index_2_base + blockIdx.y);  // the subject
		
	unsigned leader_2_number = *(leader_association_base + entry2);
		
	unsigned int position = threadIdx.x;
	
	bool *result_pointer = thrust::raw_pointer_cast( & (*(result + leader_2_number)));
	bool result_value = *result_pointer;
	
	if ( !result_value )  // die Liste muss mit 0 initialisiert sein
	{	
		// conditional compilation
		if (second_gene_entries_are_leaders)  // if the second genes entries are leaders, a different comparison strategy applies
		{
			result_value = compare_bitset_values(genes_base_1 + entry1, genes_base_2 + blockIdx.y);
		} 
		else {
			result_value = compare_bitset_values(genes_base_1 + entry1, genes_base_2 + entry2);
		}
		
		if (result_value)
		{
						
			// distribute shared memory
			char* subjstring  = s;
			char* querystring = s[blockDim.x];
			unsigned short *line_values   = (unsigned short*)&querystring[blockDim.x];
			unsigned short *intermediate  = &line_values[blockDim.x];
			bool           *non_identity  = (bool*)&intermediate[blockDim.x];
			
			// local memory
			unsigned int position  = threadIdx.x;
			unsigned int query_len = *(str_len_base + (*(index_base + entry1)));
				          subj_len = *(str_len_base + (*(index_base + entry2)));
				 
			unsigned int schwelle = epsilon_functor(query_len, subj_len);
			
			if (query_len + schwelle > subj_len || subj_len + schwelle > query_len)
			{
				if (position == 0) result_value = false;
			}
			else
			{
				char* query_string_pointer = thrust::raw_pointer_cast( & (*(str_len_base + (entry1 * string_width))));
				char* subj_string_pointer  = thrust::raw_pointer_cast( & (*(str_len_base + (entry2 *  string_width))));
				
				// copy strings into shared memory
				subjstring[position]  = subj_string_pointer[position];
				querystring[position] = query_string_pointer[position];
					
				// then calculate Hamming distance
											
				bool neq;
				if (position < query_len) neq = subjstring[position] != querystring[position]; else neq = false;
				unsigned int hamming = __syncthreads_count(neq);  // works only, if there is one block per string!!!
				
				if (onlyHamming)   // template specialization
				{
					*result_pointer = query_len == subj_len && hamming <= schwelle;
					return;
				} else 
				{
					
					if ((query_len == subj_len) && hamming <= schwelle) 
					{
						result_value = true;
					}
					else
					{							
						// prepare 1st line
						line_values[position] = position;
						
						// walk throught the characters of the query string
						for (unsigned int str2pos = 0; str2pos < query_len; ++str2pos)   
						{
						    non_identity[position] = false;	
						    intermediate[position] = str2pos + 1;
						
							if (position > 0) 
							{
								if (subjstring[position-1] == querystring[str2pos]) {
									intermediate[position] = line_values[position-1];
								}
								else
								{
							        intermediate[position] = 1 + minimum(line_values[position-1], line_values[position]);
							        non_identity[position] = true;
								}					    
							}	    
							
							__syncthreads();  // finish the above stuff and go to the next step, logically applies to all threads!!!!
							    
							unsigned int i;    
							// this code uses a partitioning of the logical strucure in the Levenshtein distance matrix
							if (position <= subj_len) 
							{
							    if (!non_identity[position]) {
									i = position + 1;
									line_values[position] = intermediate[position];
												    
								    while (non_identity[i] && i <= subj_len)
								    {
										line_values[i] = minimum(intermediate[i], line_values[i-1] + 1);
										++i;
									}
								}
							}
							__syncthreads();  //  logically applies to all threads!!!!
								
						}
						// write result
						result_value = line_values[subj_len] <= schwelle)? 1 : 0 ;
					}	
				}
			}
		}
		if (position == 0 && result_value) atomicExch(result_pointer, result_value);
	}
}


//******************************************************************************************************************************************************************************

void SequencesSet::

//******************************************************************************************************************************************************************************

void SequencesSet::determine_direct_density_reachability( thrust::device_vector<unsigned int>::iterator currentLeader, 
                                                    thrust::device_vector<unsigned int>& LeaderIndices, 
                                                    thrust::device_vector<bool>&  MaskVector)
{
	// A mask is provided and will be updated
	thrust::device_vector<unsigned int> indices1(leader_1_size);
	thrust::device_vector<unsigned int> indices2(n_entries);
	thrust::device_vector<unsigned int> 
		
	
		
	auto index_counter = thrust::make_counting_iterator(0u);
	thrust::copy_if(thrust::device, index_counter, index_counter + n_entries, 
	
	
	
	dim3 grid_dims( ... );
	
	
	
	//
	using namespace thrust;
	unsigned int stringlength = 1 + width / 32;
	stringlength *= 32;
	
	unsigned int shared_memory_usage = (sizeof(char) + sizeof(char) + sizeof(unsigned short) + sizeof(unsigned short) + sizeof(bool)) * stringlength;
	
	if (use_hamming) {
		hamming_kernel<<<n_entries, width+1, shared_memory_usage>>>		
		
		(CDR3_strings.begin(), CDR3_offsets.begin(), intermediate_result.begin(), n_entries, c_qstring_len, CDR3_threshold);
	} else 
	{
		epsilon_3_functor eps_fun(threshold);
		three_epsilon_kernel<<<n_entries, width+1, shared_memory_usage>>>(currentLeader, LeaderIndices.begin(), Strings.begin(), StringLengths.begin(), width, MaskVector.begin(), eps_fun);
	}	
}
