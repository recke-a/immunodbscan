#include "CudaDefinitions.hcu"
// #include "BitsetVector.hcu"
#include <string>
#include <vector>
#include <queue>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/fill.h>
#include <thrust/iterator/iterator_adaptor.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/remove.h>
#include <thrust/execution_policy.h>
#include <iostream>
#include <iomanip>
// #include "ImmunoDBSCAN.h"


// Constructor
// Input is the column "gene_vector" of strings with the names of VH, JH or DH genes.
// In the container, each bit of an entry corresponds to one gene name, according to "gene_reference"

template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
bitset_container::bitset_container(
	InputStringsIteratorType1 gene_vector_iter, 
	InputStringsIteratorType2 gene_reference, 
	unsigned int _n_entries, 
	unsigned int n_genes)
{
	bit_width    = (n_genes + sizeof(int)*8 - 1) / (8*sizeof(int));
	bit_wordsize = sizeof(int)*8;
	n_entries    = _n_entries;
	container    = thrust::device_vector<int>(bit_width * n_entries, 0);
		
	auto my_gene_strings = make_stringset(gene_vector_iter, _n_entries);
	unsigned int i = 0;
	for (auto ref_it = gene_reference; ref_it != gene_reference + n_genes; ++ref_it, ++i)
	{
		walk_queries(my_gene_strings.StringSet_begin(), *ref_it, my_gene_strings.get_max_strlen(), i);
	}
}

// this a container to be filled later
// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
bitset_container::bitset_container(unsigned int _n_entries, unsigned int _bit_width)
{
	bit_width    = _bit_width;
	bit_wordsize = sizeof(int)*8;
	n_entries    = _n_entries;
	container    = thrust::device_vector<int>(bit_width * _n_entries, 0);
}

//************************

unsigned int bitset_container::get_width()
{
	return bit_width;
}


unsigned int bitset_container::get_element_count()
{
	return n_entries;
}
  


// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
unsigned int bitset_container::get_offset(unsigned int xth_bit)
{
	return xth_bit / bit_wordsize;
}

// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
unsigned int bitset_container::get_bitposition(unsigned int xth_bit)
{
	return xth_bit % bit_wordsize;
}

// for communication about the container, only iterators a used
// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
typename bitset_container::entry_iterator_type bitset_container::begin_bitset()
{
	return jumping_iterator<iterator_type>(container.begin(), bit_width);
	/*auto xer = make_dual_transform_iterator(thrust::make_counting_iterator<unsigned int>(0), thrust::make_constant_iterator<unsigned int>(bit_width), thrust::multiplies<unsigned int>());
	return thrust::make_permutation_iterator(container.begin(), xer);*/
}

// to erase data when not needed
// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
void
bitset_container::erase(bitset_container::entry_iterator_type from, bitset_container::entry_iterator_type to)
{
	unsigned int distance_from = (from - begin_bitset()) * bit_width;
	unsigned int distance_to   = (to - begin_bitset()) * bit_width;
	container.erase(container.begin() + distance_from, container.begin() + distance_to);
}


// the kernel functions that deal with the bit sets are all in this file

// query is in constant memory
__device__ bool get_bit(int* wordX, unsigned int bitposition)
{
	return ((*wordX) | (1 << bitposition)) > 0;
}

__device__ void set_bit(int* wordX, unsigned int bitposition)
{
	atomicOr(wordX, 1 << bitposition);
}

// ein Kernel, um in Strings zu suchen
template <typename StringListIteratorType, typename BitsetIteratorType>
__global__
void 
find_string_occurences(
                       StringListIteratorType StringListIterator,
                       BitsetIteratorType BitsetIterator,
                       char* query_source,
                       unsigned int query_len,
                       unsigned int word_index,
                       unsigned int bit_index                   
                       )
{
	unsigned int entry_num = blockIdx.x;
	int *BitsetWordRef = thrust::raw_pointer_cast(&(*(BitsetIterator + entry_num)));
	
	extern __shared__ int s[];
	
	unsigned int alignedQuerySize = ((query_len+4-1)/4)*4;
	
	char* query_string = (char*)s;
	char* subj_string  = &query_string[alignedQuerySize];
		
	char* subj_source  = thrust::raw_pointer_cast(&thrust::get<0>(*(StringListIterator + entry_num)));
	
	for (unsigned i = threadIdx.x; i < query_len; i+= blockDim.x) query_string[i] = query_source[i];
	
	unsigned int subj_len = thrust::get<1>(*(StringListIterator + entry_num));
	for (unsigned i = threadIdx.x; i < subj_len; i+= blockDim.x) subj_string[i] = subj_source[i];
	
	__syncthreads();
	
	unsigned int position = 0;
	int hamming = 1;
	while (hamming > 0 && (position + query_len) <= subj_len)
	{
		unsigned int qpos = threadIdx.x;
		unsigned int spos = threadIdx.x + position;
		
		int neq ;
		if (qpos < query_len)
		{
			neq = subj_string[spos] != query_string[qpos];
		} 
		else 
		{
			neq = 0;
		}
		
		++position;
		hamming = __syncthreads_count(neq);
		
	}
	
	if (hamming == 0 && threadIdx.x == 0) 
	{
		set_bit(BitsetWordRef + word_index, bit_index);
		
	}
}


// This is the function that communicates with the kernel function "find_string_occurences"
// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
template <class StringSetIteratorType>
void bitset_container::walk_queries(
                       // typename StringSet<InputStringsIteratorType1>::
                       StringSetIteratorType StringListIterator,
					   const std::string& querystring, 
					   unsigned int max_strlen, 
                       unsigned int xth_bit)
{
	// first: copy querystring to constant device memory
	char* dev_query_string;
	cudaMalloc((char**)&dev_query_string, sizeof(char)*querystring.size());
	cudaMemcpy(dev_query_string, (char*)&(*querystring.begin()), sizeof(char)*querystring.size(), cudaMemcpyHostToDevice);
		
	// calculate coordinates
	unsigned int word_index = get_offset(xth_bit);
	unsigned int bit_index  = get_bitposition(xth_bit);
	
	unsigned alignedQuerySize = ((querystring.size()+4-1)/4)*4;
	unsigned alignedMaxStrlen = ((max_strlen+4-1)/4)*4;
	
	unsigned int shared_mem_usage = sizeof(char)*alignedQuerySize + sizeof(char)*alignedMaxStrlen;
	
	// then start kernel
	find_string_occurences<<<n_entries, querystring.size(), shared_mem_usage>>>(StringListIterator, begin_bitset(), dev_query_string, querystring.size(), word_index, bit_index);

    cudaFree(dev_query_string);
	
	// std::cout << "Gene " << xth_bit << "(" << word_index << "|" << bit_index << ") " << querystring << " with " << c1 << " queries and " << c2 << " positives\n"; std::cout.flush();
}


//************************************************************************************************************************
	// Version 1
	template <int sN, class first_iterator_tuple_type, class second_iterator_tuple_type, typename max_width_tuple_type> 
	__device__ 
	int compare_bitset_values(first_iterator_tuple_type bitvalues1, second_iterator_tuple_type bitvalues2, max_width_tuple_type gene_width_tuple, Int2Type<sN>)
	{
		unsigned int current_width = thrust::get<sN>(gene_width_tuple);
		
		typedef Int2Type<sN - 1> mapType;
		mapType mapInt;
		
		int cmp = 0;
		int result = compare_bitset_values(bitvalues1, bitvalues2, gene_width_tuple, mapInt);
		
		int *pointer1 = thrust::raw_pointer_cast( & (thrust::get<sN>(*bitvalues1)) );
		int *pointer2 = thrust::raw_pointer_cast( & (thrust::get<sN>(*bitvalues2)) );
		
		for (int position = threadIdx.x; position < current_width; position += blockDim.x)
		{
			cmp |= (pointer1[position] & pointer2[position]);
		}
		
		result = result && (__syncthreads_count(cmp) > 0);
		
		return result;
	}	
		
	// Version 2	
	template <class first_iterator_tuple_type, class second_iterator_tuple_type, typename max_width_tuple_type> 
	__device__ 
	int compare_bitset_values(first_iterator_tuple_type bitvalues1, second_iterator_tuple_type bitvalues2, max_width_tuple_type gene_width_tuple, Int2Type<0>)
	{
		unsigned int current_width = thrust::get<0>(gene_width_tuple);
		
		int *pointer1 = thrust::raw_pointer_cast( & (thrust::get<0>(*bitvalues1)) );
		int *pointer2 = thrust::raw_pointer_cast( & (thrust::get<0>(*bitvalues2)) );
		
		int cmp = 0;
		for (int position = threadIdx.x; position < current_width; position += blockDim.x)
		{
			cmp |= (pointer1[position] & pointer2[position]);
		}
		
		return (__syncthreads_count(cmp) > 0);
	}	
	
	
	// BaseCaller
	template <class first_iterator_tuple_type, class second_iterator_tuple_type, typename max_width_tuple_type> 
	__device__ 
	int compare_bitset_values(first_iterator_tuple_type bitvalues1, second_iterator_tuple_type bitvalues2, max_width_tuple_type gene_width_tuple)
	{
		static_assert(std::is_same<typename first_iterator_tuple_type::value_type, typename second_iterator_tuple_type::value_type>::value == true, "First and seoncd iterator value types must be equal");
		typedef Int2Type<thrust::tuple_size<typename first_iterator_tuple_type::value_type>::value-1> mapType;
		mapType map;
		
		return compare_bitset_values(bitvalues1, bitvalues2, gene_width_tuple, map);
	}
	
	

// testing function
// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
void bitset_container::example()
{
	
	std::cout << "\nExample Bitset from " << n_entries << " entries \n";
	for (unsigned int j=0; j < 10 & j < n_entries; ++j)
	{
		std::cout << j << ": ";
		for (unsigned int py = 0; py < bit_width; ++py) std::cout << std::setw(10) << std::setfill(' ') << container[bit_width*j + py];
		std::cout << std::setw(10) << std::setfill(' ') << "\nin bits = ";
		
		for (unsigned int px = 0; px < bit_width*bit_wordsize; ++px)
		{
			if (px % bit_wordsize == 0) std::cout << " ";
			char res = ( (container[bit_width*j + get_offset(px)]) & (1 << get_bitposition(px)))? '1' : '.';
			std::cout << res;
			
		}
		std::cout << "\n";
	}
}

// testing function
// template <class InputStringsIteratorType1, class InputStringsIteratorType2>	
void bitset_container::print_example_entry(unsigned int entry_index)
{
	for (unsigned int px = 0; px < bit_width * bit_wordsize; ++px)
		{
			if (px % bit_wordsize == 0) std::cout << " ";
			char res = ( (container[bit_width*entry_index + get_offset(px)]) & (1 << get_bitposition(px)))? '1' : '.';
			std::cout << res;
			
		}
		std::cout << "\n";
}



//************************************************************************************************************************

