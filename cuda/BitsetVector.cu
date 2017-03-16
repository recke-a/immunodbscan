#include "BitsetVector.hpp"
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




//************************************************************************************************************************

#define GPUERRCHK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


//************************************************************************************************************************

//__global__ void compare_bitset_entries(unsigned int query_num, bitset_vector::iterator index_iterator, 
                                  //unsigned int max_subject_num, thrust::device_vector<bool>::iterator result)
//{
	//unsigned int entry_num    = blockIdx.x * blockDim.x + threadIdx.x;
	//extern __shared__ unsigned int check[]; 
	
	//unsigned int width_pos    = threadIdx.y;
	//unsigned int subject_pos  = width_pos               + entry_num * blockDim.y;
	//unsigned int query_pos    = width_pos               + query_num *  blockDim.y;
	//unsigned int shared_pos   = width_pos               + threadIdx.x * blockDim.y;
	
	//if (entry_num < max_subject_num)
	//{	
		//check[shared_pos] = (*(index_iterator + subject_pos)) & (*(index_iterator + query_pos));   // bit-wise AND
	//}
	
	//__syncthreads();	
			
	//unsigned int n_total_threads = blockDim.y;
	    
    //while (n_total_threads > 1)
    //{
		//unsigned int half_point = (n_total_threads + 2 - 1) / 2;
		//unsigned int second = threadIdx.y + half_point;
		
		//if (second < blockDim.y)
		//{
			//check[threadIdx.y + threadIdx.x * blockDim.y] |= check[second + threadIdx.x * blockDim.y];  // bit-wise OR
		//}
		//__syncthreads();	
	
		//n_total_threads = half_point;
	//}
		
	//if (threadIdx.y == 0 & entry_num < max_subject_num)
	//{
		//*(result + entry_num) = check[threadIdx.x * blockDim.y] > 0;
	//}	
//}

//************************************************************************************************************************

bitset_vector::bitset_vector(unsigned int _n_elem, unsigned int _width) :
  wordsize(sizeof(int)*8),
  n_elem(_n_elem),
  width((_width + sizeof(int)*8 - 1) / (8*sizeof(int))),
  container( ((_width + sizeof(int)*8 - 1) / (8*sizeof(int))) * _n_elem, 0)
  {
  }
  
unsigned int bitset_vector::get_offset(unsigned int xth_bit)
{
	return xth_bit / wordsize;
}

unsigned int bitset_vector::get_bitposition(unsigned int xth_bit)
{
	return xth_bit % wordsize;
}



// thrust version


// derive repeat_iterator from iterator_adaptor
template<typename Iterator>
  class repeat_iterator
    : public thrust::iterator_adaptor< repeat_iterator<Iterator>, Iterator >
{
  public:
    typedef thrust::iterator_adaptor<
      repeat_iterator<Iterator>,
      Iterator
    > super_t;
    __host__ __device__
    repeat_iterator(const Iterator &x, int n) : super_t(x), begin(x), n(n) {}
    
    friend class thrust::iterator_core_access;
  private:
    // repeat each element of the adapted range after n steps
    unsigned int n;
    
    const Iterator begin;
    
    __host__ __device__
    typename super_t::reference dereference() const
    {
      return *(begin + (this->base() - begin) % n);
    }
};


struct key_by_width_op : public thrust::unary_function< unsigned int, unsigned int>
{
	unsigned int width;
	
	key_by_width_op(unsigned int _width) : width(_width) {}
	
	__host__ __device__ unsigned int operator()(unsigned int x)
	{
		return x / width;
	}
}; 

struct bitwise_and_to_bool_op : public thrust::unary_function< thrust::tuple<unsigned int, unsigned int>, bool>
{
	__host__ __device__ bool operator()(thrust::tuple<unsigned int, unsigned int> arg)
	{
		return arg.get<0>() & arg.get<1>();
	}
}; 

struct zipped_equal_to_op : public thrust::unary_function< thrust::tuple<unsigned int, unsigned int>, bool>
{
	__host__ __device__ bool operator()(thrust::tuple<unsigned int, unsigned int> arg)
	{
		return arg.get<0>() == arg.get<1>();
	}
}; 

	
thrust::device_vector<bool> bitset_vector::get_equivalent_entries(unsigned int query_index)
{
	thrust::device_vector<unsigned int> keys(container.size());
	thrust::device_vector<bool>	result(container.size());
	
	auto key_it   = thrust::make_transform_iterator(thrust::make_counting_iterator<unsigned int>(0), key_by_width_op(width));
	container_t::iterator query_start_iter = container.begin() + (query_index*width);
	repeat_iterator<container_t::iterator> repeat_it(query_start_iter, width); 
	auto comparing_iterator = thrust::make_zip_iterator(thrust::make_tuple(repeat_it, container.begin()));
		
	auto ende = reduce_by_key(thrust::device, 
	                          key_it, key_it + container.size(), 
						 	  thrust::make_transform_iterator(comparing_iterator, bitwise_and_to_bool_op()),
						 	  keys.begin(),
						 	  result.begin(),
						 	  thrust::equal_to<unsigned int>(),
						 	  thrust::logical_or<bool>());						 	  
	result.erase(ende.second, result.end());
	
	/*
	
	unsigned y_threads = wordsize;
	unsigned x_threads = 1024 / y_threads;
	
	unsigned n_blocks  = (n_elem+x_threads-1) / x_threads;
	unsigned dynamic_memory_size = sizeof(unsigned int) * x_threads * y_threads;
	
	dim3 blockdims(x_threads, y_threads, 1);
	
	compare_bitset_entries<<<n_blocks, blockdims, dynamic_memory_size>>>(query_index, container.begin(), n_elem, result.begin());
	*/	
		
	return result;
}	
	
	
void bitset_vector::get_equivalent_entries_with_mask(unsigned int query_index, thrust::device_vector<bool>& mask_vector)
{
	thrust::device_vector<unsigned int> keys(container.size());
	thrust::device_vector<bool>	result(container.size());
	
	auto key_it   = thrust::make_transform_iterator(thrust::make_counting_iterator<unsigned int>(0), key_by_width_op(width));
	container_t::iterator query_start_iter = container.begin() + (query_index*width);
	repeat_iterator<container_t::iterator> repeat_it(query_start_iter, width); 
	auto comparing_iterator = thrust::make_zip_iterator(thrust::make_tuple(repeat_it, container.begin()));
		
	auto ende = reduce_by_key(thrust::device, 
	                          key_it, key_it + container.size(), 
						 	  thrust::make_transform_iterator(comparing_iterator, bitwise_and_to_bool_op()),
						 	  keys.begin(),
						 	  result.begin(),
						 	  thrust::equal_to<unsigned int>(),
						 	  thrust::logical_or<bool>());						 	  
	result.erase(ende.second, result.end());
	
	thrust::transform(thrust::device, mask_vector.begin(), mask_vector.end(), result.begin(), thrust::logical_and<bool>());
}




void bitset_vector::get_identical_entries_with_mask(unsigned int query_index, thrust::device_vector<bool>& mask_vector)
{
	thrust::device_vector<unsigned int> keys(container.size());
	thrust::device_vector<bool>	result(container.size());
	
	auto key_it   = thrust::make_transform_iterator(thrust::make_counting_iterator<unsigned int>(0), key_by_width_op(width));
	container_t::iterator query_start_iter = container.begin() + (query_index*width);
	repeat_iterator<container_t::iterator> repeat_it(query_start_iter, width); 
	auto comparing_iterator = thrust::make_zip_iterator(thrust::make_tuple(repeat_it, container.begin()));
		
	auto ende = reduce_by_key(thrust::device, 
	                          key_it, key_it + container.size(), 
						 	  thrust::make_transform_iterator(comparing_iterator, zipped_equal_to_op()),
						 	  keys.begin(),
						 	  result.begin(),
						 	  thrust::equal_to<unsigned int>(),
						 	  thrust::logical_and<bool>());	 	  
	result.erase(ende.second, result.end());
	
	thrust::transform(thrust::device, mask_vector.begin(), mask_vector.end(), result.begin(), thrust::logical_and<bool>());
}


// testing function
thrust::host_vector<unsigned int> bitset_vector::get_equivalent_entry_indices(unsigned int query_index)
{
	using namespace thrust;
	device_vector<bool> equivalents = get_equivalent_entries(query_index);
	device_vector<unsigned int> indices(equivalents.size(),0);
	
	auto count_iter = make_counting_iterator(0);
	auto ende = copy_if(device, count_iter, count_iter + equivalents.size(), equivalents.begin(), indices.begin(), identity<bool>());
	host_vector<unsigned int> hv_res(indices.begin(), ende);
	return hv_res;
}


// testing function
void bitset_vector::example()
{
	
	std::cout << "\nExample Bitset from " << container.size() << " entries \n";
	for (int j=0; j < 10 & j < container.size(); ++j)
	{
		std::cout << j << ": ";
		for (int py = 0; py < width; ++py) std::cout << std::setw(10) << std::setfill(' ') << container[width*j + py];
		std::cout << std::setw(10) << std::setfill(' ') << " in bits = ";
		
		for (int px = 0; px < width*wordsize; ++px)
		{
			if (px % wordsize == 0) std::cout << " ";
			char res = ( (container[width*j + get_offset(px)]) & (1 << get_bitposition(px)))? '1' : '.';
			std::cout << res;
			
		}
		std::cout << "\n";
	}
}

// testing function
void bitset_vector::print_example_entry(unsigned int entry_index)
{
	for (int px = 0; px < width*wordsize; ++px)
		{
			if (px % wordsize == 0) std::cout << " ";
			char res = ( (container[width*entry_index + get_offset(px)]) & (1 << get_bitposition(px)))? '1' : '.';
			std::cout << res;
			
		}
}


//************************************************************************************************************************


//***** a utility kernel for finding strings and setting corresponding 

// query is in constant memory
#define MAX_QUERY_STRINGLEN 32
__constant__ char query_qstring[MAX_QUERY_STRINGLEN]; 

__global__ void find_string_and_set_bit( thrust::device_vector<char>::iterator character_iterator, 
								    thrust::device_vector<unsigned int>::iterator index_iterator, 
								    unsigned int width,
								    unsigned int offset, 
									unsigned int bitposition,
									unsigned int maxpos,
									bitset_vector::iterator bitset_vector_result)
{
	unsigned int xpos = blockIdx.x + threadIdx.x;
	
	unsigned int neq = 1;
	if (xpos < maxpos) neq = *(character_iterator+xpos) != query_qstring[threadIdx.x];
	unsigned int hamming = __syncthreads_count(neq);
	
	if (hamming == 0 && threadIdx.x == 0)   // && xpos < maxpos is not needed
	{	
		unsigned int targpos = *(index_iterator + xpos);
		atomicOr( thrust::raw_pointer_cast(&(*(bitset_vector_result + (targpos * width) + offset))), 1 << bitposition);
	}
}
                                   

//************************************************************************************************************************

StringFinderClass::StringFinderClass(std::vector<std::string>&& input, std::vector<std::string>&& _queries) :
	n_elem(input.size()),
	width(_queries.size()),		
	queries(std::move(_queries))
{
	gesamtZeichen = 0;
	for (auto str_it = input.begin(); str_it != input.end(); ++str_it) gesamtZeichen += str_it->size() + 2; 
		
	character_set = thrust::device_vector<char>(gesamtZeichen);
	index_set = thrust::device_vector<unsigned int>(gesamtZeichen);
	thrust::host_vector<char> temp_character_set(gesamtZeichen);
	thrust::host_vector<unsigned int>  temp_index_set(gesamtZeichen);
	
	auto tcs_it = temp_character_set.begin();
	auto tis_it = temp_index_set.begin();
	auto inp_it = input.begin();
	unsigned int index = 0;
	for (; inp_it != input.end(); ++inp_it, ++index, ++tcs_it, ++tis_it)
	{
		for (auto s_it = inp_it->begin(); s_it != inp_it->end(); ++s_it, ++tcs_it, ++tis_it) 
		{
			*tcs_it = *s_it;
			*tis_it = index;
		}
		*tcs_it = ' ';  // always add a delimiter to identify genes correctly!
		*tis_it = index;
		++tcs_it;
		++tis_it;
		*tcs_it = 0;
		*tis_it = index;
	}
	
	thrust::copy(temp_character_set.begin(), temp_character_set.end(), character_set.begin());
	thrust::copy(temp_index_set.begin(), temp_index_set.end(), index_set.begin());
}

bitset_vector StringFinderClass::get_bitset()
{
	bitset_vector result(n_elem, queries.size());
	
	unsigned px = 0;
	auto q_it = queries.begin();
	for (; q_it != queries.end(); ++q_it, ++px)
	{
		// transfer query to device constant memory
		char iarr[MAX_QUERY_STRINGLEN];
        thrust::fill(iarr, iarr+MAX_QUERY_STRINGLEN, 0);
		thrust::copy(q_it->begin(), q_it->end(), iarr);
		
        cudaMemcpyToSymbol(query_qstring, iarr, sizeof(char)*MAX_QUERY_STRINGLEN);
		
		find_string_and_set_bit<<<gesamtZeichen - q_it->size(), q_it->size()>>>(character_set.begin(), 
		                                                                         index_set.begin(), 
		                                                                         result.width, result.get_offset(px), result.get_bitposition(px),
		                                                                         gesamtZeichen,
		                                                                         result.container.begin());
		// GPUERRCHK(cudaPeekAtLastError());  // for debugging only
		// GPUERRCHK(cudaDeviceSynchronize());  // for debugging only
	}
	return result;
}
	
