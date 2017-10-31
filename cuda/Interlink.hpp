/*
 * InterlinkClass.hpp
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

// this file is meant to perform the communication between C++ and CUDA code

#pragma once
#include <string>
#include <memory>
#include <thrust/tuple.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/permutation_iterator.h>

/*
#ifdef __CUDACC__
	#include "BitsetVector.hpp"
	#include "StringOperations.hpp"
#endif	
*/

typedef thrust::host_vector<std::string>::iterator NormalStringVectorIteratorType;
typedef thrust::tuple<NormalStringVectorIteratorType,NormalStringVectorIteratorType,NormalStringVectorIteratorType> NormalStringVectorIteratorTupleType;
typedef thrust::zip_iterator<NormalStringVectorIteratorTupleType> NormalStringVectorZipIteratorType;


class InterlinkClass
{	
	thrust::host_vector<unsigned int> clusters, uniques, copies;
	thrust::host_vector<int> cores, vbitset;
	unsigned int n_uniques;		
	unsigned int v_bitset_width;
			
	public:
	InterlinkClass(NormalStringVectorZipIteratorType InputStrings, unsigned int _n_entries, 
			  NormalStringVectorIteratorType V_genes_reference, unsigned int V_gene_variant_count,
			  NormalStringVectorIteratorType J_genes_reference, unsigned int J_gene_variant_count,
			  float _threshold, unsigned int _minPts) ;
	
	thrust::host_vector<unsigned int>::iterator get_clusters_iterator();
	thrust::host_vector<int>::iterator get_cores_iterator();
	thrust::host_vector<unsigned int>::iterator get_uniques_iterator();
	thrust::host_vector<unsigned int>::iterator get_copies_iterator();
	unsigned int get_n_uniques();
	
	thrust::host_vector<int>::iterator v_gene_container_begin();
	thrust::host_vector<int>::iterator v_gene_container_end();
	unsigned int get_v_bitset_width() { return v_bitset_width; }
	
};
