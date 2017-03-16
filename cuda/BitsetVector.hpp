
#pragma once
#include <string>
#include <vector>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/iterator_adaptor.h>

//************************************************************************************************************************

struct bitset_vector
{
	typedef thrust::device_vector<int> container_t;
	typedef container_t::iterator iterator;
	container_t container;
    unsigned int n_elem;
    unsigned int width;
    unsigned int wordsize;
    
    bitset_vector(unsigned int _n_elem, unsigned int _width);
      
    unsigned int get_offset(unsigned int xth_bit);
	
	unsigned int get_bitposition(unsigned int xth_bit);		
    
    thrust::device_vector<bool> get_equivalent_entries(unsigned int query_index);
    
	// testing functions
    thrust::host_vector<unsigned int> get_equivalent_entry_indices(unsigned int query_index);
    void get_equivalent_entries_with_mask(unsigned int query_index, thrust::device_vector<bool>& mask_vector);
    void get_identical_entries_with_mask(unsigned int query_index, thrust::device_vector<bool>& mask_vector);
    
    void example();
    void print_example_entry(unsigned int entry_index);
      
};  
	

//************************************************************************************************************************

// derive jumping iterator from iterator_adaptor
template<typename Iterator>
  class jumping_iterator
    : public thrust::iterator_adaptor< jumping_iterator<Iterator>, Iterator >
{
  public:
    typedef thrust::iterator_adaptor<
      jumping_iterator<Iterator>,
      Iterator
    > super_t;
    
    __host__ __device__
    jumping_iterator(const Iterator &x, int n) : super_t(x), begin(x), n(n) {}
    
    friend class thrust::iterator_core_access;
    
    __host__ __device__
    unsigned int get_jumping_length() const
    {
		return n;
	}
    
  private:
    // repeat each element of the adapted range after n steps
    unsigned int n;
    
    const Iterator begin;
    
    __host__ __device__
    typename super_t::reference dereference() const
    {
      return *(begin + (this->base() - begin) * n);
    }
};

//************************************************************************************************************************


class StringFinderClass {

	thrust::device_vector<char> character_set;
	thrust::device_vector<unsigned int> index_set;
	std::vector<std::string> queries;
	unsigned int n_elem;
	unsigned int width;
	unsigned int gesamtZeichen;
	
    public:
    // This function walks through a Sequencing_Container to identify the genes

	StringFinderClass(std::vector<std::string>&& input, std::vector<std::string>&& _queries);
	
	bitset_vector get_bitset();
};

