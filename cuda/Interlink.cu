/*
 * InterlinkClass.cu
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


#include "Interlink.hpp" 
#include "BitsetVector.hcu"
#include "DataContainer.hcu"
#include "PSDBSCAN.hcu"
#include <thrust/iterator/constant_iterator.h>
#include "StringOperations.hcu"
#include <thrust/iterator/counting_iterator.h>
 
 	
InterlinkClass::InterlinkClass(NormalStringVectorZipIteratorType InputStrings, unsigned int _n_entries, 
  NormalStringVectorIteratorType V_genes_reference, unsigned int V_gene_variant_count,
  NormalStringVectorIteratorType J_genes_reference, unsigned int J_gene_variant_count,
  float _threshold, unsigned int _minPts, bool use_only_hamming, bool eval_trees) 
  :
  clusters(_n_entries, 0),
  uniques(_n_entries, 0), 
  copies (_n_entries, 0),
  cores(_n_entries, 0)
{
	bitset_container V_bitset(thrust::get<0>(InputStrings.get_iterator_tuple()), V_genes_reference, _n_entries, V_gene_variant_count);
	vbitset.assign(V_bitset.begin_container(), V_bitset.end_container());
	v_bitset_width = V_bitset.get_width();
	
	// V_bitset.example();
	
	bitset_container J_bitset(thrust::get<1>(InputStrings.get_iterator_tuple()), J_genes_reference, _n_entries, J_gene_variant_count);
		
	// J_bitset.example();	
		
	auto CDR3strings = std::move(make_stringset(thrust::get<2>(InputStrings.get_iterator_tuple()), _n_entries));
		
	auto copies_iterator = thrust::make_constant_iterator<unsigned int>(0u);
		
	using bitset1_iter = decltype(V_bitset.begin_bitset());
	using bitset2_iter = decltype(J_bitset.begin_bitset());
	using string_iter  = decltype(CDR3strings.StringSet_begin());
	using copies_iter  = decltype(copies_iterator);
	DataManagementClass<thrust::zip_iterator<thrust::tuple<bitset1_iter, bitset2_iter> >, thrust::tuple<unsigned int, unsigned int>,
		string_iter, copies_iter> 		
			my_container(thrust::make_zip_iterator(thrust::make_tuple(V_bitset.begin_bitset(), J_bitset.begin_bitset())),
						thrust::make_tuple(V_bitset.get_width(), J_bitset.get_width()),
						CDR3strings.StringSet_begin(),
						copies_iterator,
						_n_entries,
						CDR3strings.get_max_strlen()
						);
	
		
	using dbscan_container_type = decltype(my_container);	

	DBSCAN<dbscan_container_type> runner(my_container, _threshold, _minPts);
	if  (use_only_hamming) runner.use_only_hamming();
	if (eval_trees) runner.set_wish_for_tree_evaluation();
	runner.doScan();
	
	 
	thrust::copy(runner.get_cluster_ids(), runner.get_cluster_ids() + _n_entries, clusters.begin());
	
	
	thrust::copy(runner.get_cores(), runner.get_cores() + _n_entries, cores.begin());
	
	
	n_uniques = runner.get_n_uniques();
	
	
	thrust::copy(runner.get_uniques_iterator(), runner.get_uniques_iterator() + n_uniques, uniques.begin());	
	
	
	thrust::copy(runner.get_copies_iterator(), runner.get_copies_iterator() + _n_entries, copies.begin());
	
}

	
thrust::host_vector<unsigned int>::iterator InterlinkClass::get_clusters_iterator()
{
	return clusters.begin();
}

thrust::host_vector<int>::iterator InterlinkClass::get_cores_iterator()
{
	return cores.begin();
}

thrust::host_vector<unsigned int>::iterator InterlinkClass::get_uniques_iterator()
{
	return uniques.begin();
}

thrust::host_vector<unsigned int>::iterator InterlinkClass::get_copies_iterator()
{
	return copies.begin();
}


thrust::host_vector<int>::iterator InterlinkClass::v_gene_container_begin()
{
	return vbitset.begin();
}

thrust::host_vector<int>::iterator InterlinkClass::v_gene_container_end()
{
	return vbitset.end();
}


unsigned int InterlinkClass::get_n_uniques()
{
	return n_uniques;
}

