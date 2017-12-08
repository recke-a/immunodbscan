/*
 * csvreader.hpp
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
#include <stdexcept>
#include <string>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <vector>
#include <algorithm>
#include <future>
#include <memory>


// this class reads csv data and stores it into a TableHolder

class CSVReaderClass
{
	//*************************************************************************************
	const char separator;
	const char endline;
	const char quotation;
	
	std::ifstream file;
	std::string filename;	
	
	
	std::unique_ptr<std::vector<int> > column_selector_order_pt;
	
	//*************************************************************************************	
	
	template <typename col_iterator_type>
	void process_line(std::string&& dataline, unsigned int ref_line, col_iterator_type col_it)
	{
		boost::tokenizer<boost::escaped_list_separator<char> > tk(dataline, boost::escaped_list_separator<char>('\\', separator, quotation));
		(*col_it)[ref_line] = filename;
		
		auto scatter_it = column_selector_order_pt->begin();
	    auto tok_it     = tk.begin();
	    unsigned int n = 0;
	    for (; tok_it != tk.end(); ++tok_it, ++scatter_it, ++n)
	    {
			if (*scatter_it != -1) 
			{
				(*(col_it + (*scatter_it)))[ref_line] = *tok_it;
			}
		}  
	}
			
	//*************************************************************************************	
	
	void purify_line(std::string& line)
	{
		if ( line.size() && line[line.size()-1] == '\r' ) 
		{
           line = line.substr( 0, line.size() - 1 );
       }
	}
	
	//*************************************************************************************
			
	public:
	
	CSVReaderClass(const std::string& _filename, const char _separator, const char _endline, const char _quotation) :
		separator(_separator), 
	    endline(_endline),
	    quotation(_quotation),
	    file(_filename),
	    filename(_filename)
	{
		if(!file.is_open())
		{
		  throw std::runtime_error("Cannot read from file " + filename);
		}
	}
	
	
	//*************************************************************************************
			
	template <typename header_strings_iterator_type>
	void sync_headerline(header_strings_iterator_type header_strings_iterator_begin, header_strings_iterator_type header_strings_iterator_end)
	{
		std::string line;
        std::getline(file, line, file.widen(endline));
        purify_line(line);
        
        boost::tokenizer<boost::escaped_list_separator<char> > tk(line, boost::escaped_list_separator<char>('\\', separator, quotation));
        
        // transfer tokened strings to better container
        std::vector<std::string> headerline(tk.begin(), tk.end());
        
        // set order vector to -1, with -1 = n.a.
        column_selector_order_pt.reset(new std::vector<int>(headerline.size(),-1));
        
        // A bit messy, but here I set a first column to be the filename 
        int i = 1;
        auto csl_it_begin = column_selector_order_pt.get()->begin();
        for (auto it = header_strings_iterator_begin+1; it != header_strings_iterator_end; ++it, ++i)
        {
			// compare list of selected headers with headers in file
			auto pos = std::find(headerline.begin(), headerline.end(), *it); 
			
			// exit with noise, if the selected header is not found
			if (pos == headerline.end()) {
				std::cout << "\nError while checking table headers. Column \"" << (*it) << "\" was not found. These headers are available: \n***\n";
				std::for_each(headerline.begin(), headerline.end(), [&](const std::string& hl) { std::cout << "Länge = " << hl.size() << " \"" << hl << "\"\n"; });
				std::cout << "***\n"; std::cout.flush();
				throw std::runtime_error("Required column " + (*it) + " not found in " + filename);
			}
			
			// change into an index
			int ind_pos = pos - headerline.begin();
			
			// and save it for later
			*(csl_it_begin + ind_pos) = i;
		}
	}	
		
	//*************************************************************************************	
	
	template <typename col_iterator_type>
	void read_full_file(col_iterator_type col_iterator_begin, col_iterator_type col_iterator_end)
	{
		const unsigned int BLOCK_SIZE = 10000;
	
		std::string line;
		
		unsigned int cur_n_rows  = (*col_iterator_begin).size();
		unsigned int line_number = cur_n_rows;
				
		while (!file.eof()) {
			
			if (line_number == cur_n_rows)
			{
				std::for_each(col_iterator_begin, col_iterator_end, [&](typename col_iterator_type::value_type &xcol){ xcol.resize( cur_n_rows + BLOCK_SIZE); });
				cur_n_rows += BLOCK_SIZE;
			}
		
			std::vector<std::future<void> > tasklist;
			
			for (unsigned int fut_i = 0; fut_i < BLOCK_SIZE && !file.eof(); ++fut_i, ++line_number)
			{
				std::getline(file, line, file.widen(endline));
				purify_line(line);
				
				// spawn tasks
				tasklist.push_back( std::async( std::launch::async,
					[&](std::string&& data, unsigned int ref_line, col_iterator_type col_it) 
					{ 
						return process_line(std::forward<std::string>(data), ref_line, col_it);
					},
					std::move(line),
					line_number,
					col_iterator_begin )
				);
				
				// wait for all results to become available	
				std::for_each(tasklist.begin(), tasklist.end(), [&](std::future<void> &single_task) { single_task.wait(); });
			}		
			
		}
		
		// resize vectors to actual size
		if (line_number < cur_n_rows)
		{
			std::for_each(col_iterator_begin, col_iterator_end, [&](typename col_iterator_type::value_type &xcol){ xcol.resize(line_number); });
		}
	}
};
