/*
 * csvwriter.hpp
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
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include <fstream>
 
 class CSVWriterClass
 {
    //*************************************************************************************
	
	public:
	
	template <typename input_iterator_type_1, typename input_iterator_type_2>
	CSVWriterClass(const input_iterator_type_1 table_input_begin, 
				   const input_iterator_type_1 table_input_end, 
				   const input_iterator_type_2 colnames_begin,
				   const input_iterator_type_2 colnames_end,
				   const std::string& filename, 
				   const char _separator, const char _endline, const char _quotation) 
	{	
		std::ofstream out(filename, std::ofstream::out);
		// header
		auto ith = colnames_begin;
		while (ith != colnames_end) 
		{
			out << (*ith);
			++ith;
			if (ith != colnames_end) out << _separator;
		}
		out << _endline;
		
		// contents
		for (unsigned row = 0; row < (*table_input_begin).size(); ++row)
		{
			auto it = table_input_begin; 
			while (it != table_input_end)
			{
				if (_quotation != '\0')
					out << _quotation << (*it).at(row) << _quotation;
				else
					out << (*it).at(row);
					
				++it;
				if (it != table_input_end) out << _separator;
			}
			out << _endline;
		}
		out.close();
	}
	 
 };
 
