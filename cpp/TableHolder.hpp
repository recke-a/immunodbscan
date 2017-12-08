/*
 * unbenannt
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

#pragma once
#include <string>
#include <iostream>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <boost/iterator/permutation_iterator.hpp>

#include "ConfigHolder.hpp"
#include "csvreader.hpp"
#include "csvwriter.hpp"
#include "BenchmarkSingleton.hpp"
#include "../cuda/Interlink.hpp"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>

// This class implements the container for the data set

class TableHolderClass
{	
	ConfigHolderClass& ConfigHolder;
	std::vector<std::string> Col_Names;
	std::vector<std::vector<std::string> > StringTableContainer;
	
	char sep, endline, quotation;
			
	std::string single_filename;
		
	template <typename outtype>
	typename std::enable_if<!std::is_same<outtype, std::string>::value, std::string>::type make_output_string(outtype x) const
	{
		return std::to_string(x);
	}	
	
	template <typename outtype>
	typename std::enable_if<std::is_same<outtype, std::string>::value, std::string>::type make_output_string(outtype x) const
	{
		return x;
	}
		
	public:
	
	//*************************************************************************************
	
	TableHolderClass(std::vector<std::string> input_files, ConfigHolderClass& _ConfigHolder) :
		ConfigHolder(_ConfigHolder),
		Col_Names(_ConfigHolder.get_input_column_count() + 1),
		StringTableContainer(_ConfigHolder.get_input_column_count() + 1)
	{
		sep = ConfigHolder.getSpecialCharacter(ConfigHolderClass::fsInput, ConfigHolderClass::separator);
		endline = ConfigHolder.getSpecialCharacter(ConfigHolderClass::fsInput, ConfigHolderClass::endline);
		quotation = ConfigHolder.getSpecialCharacter(ConfigHolderClass::fsInput, ConfigHolderClass::quotation);
		
		Col_Names[0] = "file_origin";
		ConfigHolder.get_input_header_denominators(Col_Names.begin()+1);
		
		single_filename = *(input_files.begin());
		
		// now ... read the input files (plural, if merging is wanted)
		std::for_each(input_files.begin(), input_files.end(), [&](const std::string& filename) { read_input_file(filename); });
	}
		
	
	//*************************************************************************************
		
	void read_input_file(const std::string& input_file)
	{
		// start reader
		CSVReaderClass CSVReader(input_file, sep, endline, quotation);
				
		// encapsulate this within this function
		std::vector<std::string> CSV_Col_Names(ConfigHolder.get_input_column_count()+1,"");  // the first is left empty!
		ConfigHolder.get_input_header_strings(CSV_Col_Names.begin()+1);
		
		// first, check headerline for problems and sync for full-speed reading
		CSVReader.sync_headerline(CSV_Col_Names.begin(), CSV_Col_Names.end());
		
		// then go for reading
		CSVReader.read_full_file(StringTableContainer.begin(), StringTableContainer.end());
	}
		
		
	//*************************************************************************************	
		
	unsigned int get_total_rows_count()
	{
		return StringTableContainer[0].size();
	}
	
	//*************************************************************************************	
	
	unsigned int get_total_col_count()
	{
		return StringTableContainer.size();
	}
	
	//*************************************************************************************	
	
	template <typename input_iterator_type>
	void append_column(const std::string& colname, input_iterator_type input_iterator_begin, input_iterator_type input_iterator_end)
	{
		std::vector<std::string> app_col(get_total_rows_count(), "");
		std::transform(input_iterator_begin, input_iterator_end, app_col.begin(), [&](typename input_iterator_type::value_type x) { return make_output_string(x); });
		Col_Names.push_back(colname);
		StringTableContainer.emplace_back(std::move(app_col));
	}
	
	//*************************************************************************************	
	
	std::vector<std::string>::iterator get_column_iterator(const std::string& colname)
	{
		auto whichcol_it = std::find(Col_Names.begin(), Col_Names.end(), colname);
		if (whichcol_it == Col_Names.end()) throw std::runtime_error("Column " + colname + " not contained in TableHolder!"); 
		unsigned int pos = whichcol_it - Col_Names.begin();
		return StringTableContainer[pos].begin();
	}
	
	//*************************************************************************************	
	
	template <typename label_iterator_type>
	void delete_rows_by_label(label_iterator_type begin_label, label_iterator_type end_label)
	{
		unsigned int n_deleted = std::accumulate(begin_label, end_label, 0);
		unsigned int new_total_rows = get_total_rows_count() - n_deleted;
		
		BenchmarkSingleton::Instance().Message("Deleting " + std::to_string(n_deleted) + " rows");
			
		std::for_each(StringTableContainer.begin(), StringTableContainer.end(), [&](std::vector<std::string>& stvec)
			{
				auto it = begin_label;
				auto sit = stvec.begin();
				auto tit = stvec.begin();
				
				for (; it != end_label; ++it, ++sit)
				{
					*tit = *sit;
					if (! (*it)) ++tit;
				}
				
				stvec.resize(new_total_rows);
			});		
	}
	
	//*************************************************************************************	
		
	void delete_empty_values(const std::string& colname)
	{
		std::vector<int> to_be_deleted(get_total_rows_count(), 0);
		
		// searches for "" to mark row as deletable
		auto whichcol_it = std::find(Col_Names.begin(), Col_Names.end(), colname);
		if (whichcol_it == Col_Names.end()) throw std::runtime_error("Column " + colname + " not contained in TableHolder!"); 
		unsigned int pos = whichcol_it - Col_Names.begin();
		
		std::transform(StringTableContainer[pos].begin(), 
					   StringTableContainer[pos].end(),
					   to_be_deleted.begin(),
					   [&](const std::string& st) {
						   return st.compare("") == 0;
					   });
	
		BenchmarkSingleton::Instance().Message("Deleting empty rows for column " + colname);
	
		delete_rows_by_label(to_be_deleted.begin(), to_be_deleted.end());
	}	
	
	//*************************************************************************************	
	
	void delete_noise_rows()
	{
		std::vector<int> to_be_deleted(get_total_rows_count(), 0);
		
		// searches for "" to mark row as deletable
		auto whichcol_it = std::find(Col_Names.begin(), Col_Names.end(), "clone_id");
		unsigned int pos = whichcol_it - Col_Names.begin();
		
		std::transform(StringTableContainer[pos].begin(), 
					   StringTableContainer[pos].end(),
					   to_be_deleted.begin(),
					   [&](const std::string& st) {
						   return st.compare("0") == 0;
					   });
	
		BenchmarkSingleton::Instance().Message("Deleting noise points");
		delete_rows_by_label(to_be_deleted.begin(), to_be_deleted.end());
	}	
	
	
	//*************************************************************************************	
	
	void write_to_file()
	{
		char wsep = ConfigHolder.getSpecialCharacter(ConfigHolderClass::fsOutput, ConfigHolderClass::separator);
		char wendline = ConfigHolder.getSpecialCharacter(ConfigHolderClass::fsOutput, ConfigHolderClass::endline);
		char wquotation = ConfigHolder.getSpecialCharacter(ConfigHolderClass::fsOutput, ConfigHolderClass::quotation);
		
		// do we need to remove noise?
		if (ConfigHolder.remove_noise_from_output())
		{
			delete_noise_rows();
		}
		
		// the output columns that are not available are added as empty columns
		
		unsigned int n_output_columns = ConfigHolder.get_output_column_count();
		
		// read headernames
		std::vector<std::string> output_headers(n_output_columns);
		ConfigHolder.get_output_header_strings(output_headers.begin());
		
		// read colname denominators
		std::vector<std::string> selected_table_cols(n_output_columns);
		ConfigHolder.get_output_header_denominators(selected_table_cols.begin());
		
		// compare selected_table_cols with data in this container
		std::vector<unsigned int> new_order_of_columns(n_output_columns);
		
		std::for_each(selected_table_cols.begin(), selected_table_cols.end(), 
			[&](const std::string& st) 
			{
				if (std::find(Col_Names.begin(), Col_Names.end(), st) == Col_Names.end())
				{
					std::vector<std::string> emptyvec(get_total_rows_count(), "");
					append_column(st, emptyvec.begin(), emptyvec.end());
				}
			});
		
		// copy refs of StringTable contents to target vector
		std::transform(selected_table_cols.begin(), selected_table_cols.end(), new_order_of_columns.begin(),
		    [&](const std::string& st) 
			{
				auto t_it = std::find(Col_Names.begin(), Col_Names.end(), st);
				unsigned int t_pos = t_it - Col_Names.begin();
				return t_pos;
			});
			
		auto new_order_iterator_begin = boost::make_permutation_iterator( StringTableContainer.begin(), new_order_of_columns.begin() );
		
		BenchmarkSingleton::Instance().Message("Writing to file " + ConfigHolder.create_filename(single_filename, ConfigHolderClass::fsOutput));
		
		CSVWriterClass(new_order_iterator_begin, new_order_iterator_begin + n_output_columns, 
					   output_headers.begin(), output_headers.end(),
					   ConfigHolder.create_filename(single_filename, ConfigHolderClass::fsOutput),
					   wsep, wendline, wquotation);
	}	
	
	//**************************************************************************************************
	
	void do_the_clustering()
	{
		BenchmarkSingleton::Instance().StartTimer("Deleting Empty Rows");
	
		std::string first_gene_colname = ConfigHolder.get_DBSCAN_column_name(ConfigHolderClass::dbFirstGene);
		std::string second_gene_colname = ConfigHolder.get_DBSCAN_column_name(ConfigHolderClass::dbSecondGene);
		std::string sequence_colname = ConfigHolder.get_DBSCAN_column_name(ConfigHolderClass::dbSequence);
	
		if (sequence_colname.compare("") != 0)   delete_empty_values(sequence_colname);
		if (first_gene_colname.compare("") != 0) delete_empty_values(first_gene_colname);
		if (second_gene_colname.compare("") != 0) delete_empty_values(second_gene_colname);
			
		BenchmarkSingleton::Instance().StopTimer("Deleting Empty Rows");
		BenchmarkSingleton::Instance().StartTimer("Transferring data to device");		
		// first gene
		
		thrust::host_vector<std::string> v_entry_string, v_gene_names;
		
		if (first_gene_colname.compare("") == 0) {
			v_entry_string.assign(get_total_rows_count(), "V");  // otherwise, I get strange warnings
			v_gene_names.assign(1, "V");
		}
		else
		{
			v_entry_string.resize(get_total_rows_count());
			thrust::copy(thrust::host, get_column_iterator(first_gene_colname), get_column_iterator(first_gene_colname) + get_total_rows_count(), v_entry_string.begin());
   			Gene_Reader VGeneReader(ConfigHolder.get_DBSCAN_column_ref_file(ConfigHolderClass::dbFirstGene));
			v_gene_names.assign(VGeneReader.begin(), VGeneReader.end());
    	}
		
		// second gene
		
		thrust::host_vector<std::string> j_entry_string, j_gene_names;
		
		std::cout.flush();
		
		if (second_gene_colname.compare("") == 0) {
			j_entry_string.assign(get_total_rows_count(), "V");  // otherwise, I get strange warnings
			j_gene_names.assign(1, "V");
		}
		else
		{
			j_entry_string.resize(get_total_rows_count());
			thrust::copy(thrust::host, get_column_iterator(second_gene_colname), get_column_iterator(second_gene_colname)+get_total_rows_count(), j_entry_string.begin());
   			Gene_Reader JGeneReader(ConfigHolder.get_DBSCAN_column_ref_file(ConfigHolderClass::dbSecondGene));
			j_gene_names.assign(JGeneReader.begin(), JGeneReader.end());
		}
		
		// sequence
		thrust::host_vector<std::string> sequence_strings;
	
		if (sequence_colname.compare("") == 0) {
			sequence_strings.assign(get_total_rows_count(), "Gene");  // otherwise, I get strange warnings
		}
		else
		{
			sequence_strings.assign(get_column_iterator(sequence_colname), get_column_iterator(sequence_colname)+get_total_rows_count());
		}
		
		BenchmarkSingleton::Instance().StopTimer("Transferring data to device");		
	
		// ready for more
		NormalStringVectorZipIteratorType 
		  InputIterators = thrust::make_zip_iterator(thrust::make_tuple(v_entry_string.begin(), j_entry_string.begin(), sequence_strings.begin()));
		
		InterlinkClass Interlink(InputIterators, get_total_rows_count(), 
			v_gene_names.begin(), v_gene_names.size(),
			j_gene_names.begin(), j_gene_names.size(),
			ConfigHolder.get_DBSCAN_threshold(), ConfigHolder.get_DBSCAN_minPts(), ConfigHolder.is_set_only_hamming(), ConfigHolder.get_tree_depth_statistics());
		
		append_column("clone_id", Interlink.get_clusters_iterator(), Interlink.get_clusters_iterator() + get_total_rows_count());
		append_column("copies", Interlink.get_copies_iterator(), Interlink.get_copies_iterator() + get_total_rows_count());
		append_column("cores", Interlink.get_cores_iterator(), Interlink.get_cores_iterator() + get_total_rows_count() );
	}

	
	
};
