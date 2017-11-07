/*
 * ConfigHolder.hpp
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <vector>

namespace pt = boost::property_tree;  

class ConfigHolderClass
{
	pt::ptree tree;
	
	//*************************************************************************************************
	
	std::string create_filename(const std::string& pathname,const std::string& prefix, const std::string& mainname, const std::string& suffix, const std::string& filetype)
	{
		return pathname + prefix + mainname + suffix + filetype;
	}
	
	
	public:
	//*************************************************************************************************
	
	enum FileSelectorType { fsInput, fsOutput, fsBenchmark };
	enum SpecialCharacterType { separator, decimal, endline, quotation };
	enum DBSCANColumnSelectionType { dbFirstGene, dbSecondGene, dbSequence };
	
	//*************************************************************************************************
	
	ConfigHolderClass(const std::string& config_file)
	{
		pt::read_info(config_file, tree);
	}
	
	
	//*************************************************************************************************
	
	char getSpecialCharacter(FileSelectorType FileSelector, SpecialCharacterType SpecialCharacter)
	{
		std::string file_which, sc_which;
		switch (FileSelector)
		{
			case fsInput:     file_which = "input_settings";  break;
			case fsOutput:    file_which = "output_settings"; break;
			case fsBenchmark: file_which = "benchmark_output_settings"; break;
		};
		
		switch (SpecialCharacter)
		{
			case separator:	sc_which = "separator"; break;
			case decimal:	sc_which = "decimal";   break;
			case endline:	sc_which = "endline";   break;
			case quotation: sc_which = "quotation"; break;
		};
		
		
		return tree.get_child(file_which).get<char>(sc_which);
	}
	
	//*************************************************************************************************	
	
	std::string create_filename(const std::string& mainname, FileSelectorType FileSelector)
	{
		std::string pathname, prefix, suffix, filetype, strippedmainname;
		
		pathname = mainname.substr(0, mainname.find_last_of("/")+1);
		strippedmainname = mainname.substr(0, mainname.find_last_of("."));
		strippedmainname = strippedmainname.substr(strippedmainname.find_last_of("/")+1, strippedmainname.length());
				
		switch (FileSelector)
		{
			case fsOutput: 
				prefix   = tree.get<std::string>("output_prefix");
				suffix   = tree.get<std::string>("output_suffix");
				filetype = tree.get<std::string>("output_filetype");
				break;
			case fsBenchmark: 
				prefix   = tree.get<std::string>("benchmark_output_prefix");
				suffix   = tree.get<std::string>("benchmark_output_suffix");
				filetype = tree.get<std::string>("benchmark_output_filetype");
				break;
			default: throw std::runtime_error("Bad File Selection!\n");
		}
		
		return create_filename(pathname, prefix, strippedmainname, suffix, filetype);
	}
	
	//*************************************************************************************************	
	
	unsigned int get_input_column_count()
	{
		auto parent_node = tree.get_child("select_columns_to_read");		
		return parent_node.size();
	}
		
	template < typename OutputIteratorClass >	
	void get_input_header_strings(OutputIteratorClass OutputIterator)
	{
		auto parent_node = tree.get_child("select_columns_to_read");
		using outtype = typename OutputIteratorClass::value_type;
		static_assert(std::is_same< outtype, std::string >::value, "Type mismatch in input header vector!\n");
		
		using iterator_value_type = decltype( *(parent_node.begin()) );
		
		std::transform(parent_node.begin(), parent_node.end(), OutputIterator, [&](iterator_value_type x) { return x.second.data(); });		
	} 
	
	template < typename OutputIteratorClass >	
	void get_input_header_denominators(OutputIteratorClass OutputIterator)
	{
		auto parent_node = tree.get_child("select_columns_to_read");
		using outtype = typename OutputIteratorClass::value_type;
		static_assert(std::is_same< outtype, std::string >::value, "Type mismatch in input header vector!\n");
		
		using iterator_value_type = decltype( *(parent_node.begin()) );
		
		std::transform(parent_node.begin(), parent_node.end(), OutputIterator, [&](iterator_value_type x) { return x.first; });		
	}
	
	//*************************************************************************************************	
	
	unsigned int get_output_column_count()
	{
		auto parent_node = tree.get_child("select_columns_to_write");		
		return parent_node.size();
	}
		
	template < typename OutputIteratorClass >	
	void get_output_header_strings(OutputIteratorClass OutputIterator)
	{
		auto parent_node = tree.get_child("select_columns_to_write");
		using outtype = typename OutputIteratorClass::value_type;
		static_assert(std::is_same< outtype, std::string >::value, "Type mismatch in input header vector!\n");
		
		using iterator_value_type = decltype( *(parent_node.begin()) );
		
		std::transform(parent_node.begin(), parent_node.end(), OutputIterator, [&](iterator_value_type x) { return x.second.data(); });		
	} 
	
	template < typename OutputIteratorClass >	
	void get_output_header_denominators(OutputIteratorClass OutputIterator)
	{
		auto parent_node = tree.get_child("select_columns_to_write");
		using outtype = typename OutputIteratorClass::value_type;
		static_assert(std::is_same< outtype, std::string >::value, "Type mismatch in input header vector!\n");
		
		using iterator_value_type = decltype( *(parent_node.begin()) );
		
		std::transform(parent_node.begin(), parent_node.end(), OutputIterator, [&](iterator_value_type x) { return x.first; });		
	}
	
	//*************************************************************************************************	
	
	float get_DBSCAN_threshold()
	{
		return tree.get<float>("maximum_relative_edit_distance");
	}
	
	unsigned int  get_DBSCAN_minPts()
	{
		return tree.get<unsigned int>("minimum_number_of_neighbors");
	}

	bool is_set_only_hamming()
	{
		return tree.get("use_hamming", false);
	}

	bool get_tree_depth_statistics()
	{
		return tree.get("get_tree_depth_statistics",false);
	}

	//*************************************************************************************************	
	
	bool remove_noise_from_output()
	{
		return tree.get<bool>("remove_noise_from_output");
	}

	//*************************************************************************************************	
	
	std::string get_DBSCAN_column_name(DBSCANColumnSelectionType what)
	{
		std::string result;
		
		try {
			switch (what)
			{
				case dbFirstGene: 
					result = tree.get_child("gene_columns").get_child("first_gene").get<std::string>("name"); break;
			
				case dbSecondGene: 
					result = tree.get_child("gene_columns").get_child("second_gene").get<std::string>("name"); break;
			
				case dbSequence:
					result = tree.get_child("sequence_columns").get<std::string>("name"); break;
			}
		} 
		catch(const pt::ptree_bad_path &e)
		{		
			result = "";
		}
		
		return result;
	}
		
		
	std::string get_DBSCAN_column_ref_file(DBSCANColumnSelectionType what)
	{
		std::string result;
		
		try {
			switch (what)
			{
				case dbFirstGene: 
					result = tree.get_child("gene_columns").get_child("first_gene").get("ref",""); break;
			
				case dbSecondGene: 
					result = tree.get_child("gene_columns").get_child("second_gene").get("ref",""); break;
					
				default:  
					throw std::runtime_error("Sequence does not have a reference file");
			}
		}
		
		catch(const pt::ptree_bad_path &e)
		{		
			result = "";
		}
		return result;
	}
};
