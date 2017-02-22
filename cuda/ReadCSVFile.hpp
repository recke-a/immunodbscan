/*
 * ReadCSVFile.hpp
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

#include <string>
#include <fstream>
#include "boost/tokenizer.hpp"
#include <vector>
#include <algorithm>


class Sequencing_Container
{
	std::vector<std::vector <std::string>> file_table;
	std::vector<std::string> headerline;
	
	std::vector<unsigned int> selected_cols;   // order: SEQUENCE_NUMBER, SEQUENCE_ID, V-GENE, D-GENE, J-GENE, CDR3
	std::vector<int> scatter_cols;   // order: SEQUENCE_NUMBER, SEQUENCE_ID, V-GENE, D-GENE, J-GENE, CDR3
	
	std::ifstream       file;
	char colsep;	
		
	void get_headerline()
	{
		std::string line;
        std::getline(file, line);
        
        boost::tokenizer<boost::escaped_list_separator<char> > tk(line, boost::escaped_list_separator<char>('\\', colsep, '\"'));
                
        headerline.assign(tk.begin(), tk.end());
	}
	
	// This function creates the map scatter_cols from selected_cols
	void invert_scatter()
	{
		scatter_cols.assign(headerline.size(), -1);
		for (int i = 0; i < selected_cols.size(); ++i) scatter_cols[selected_cols[i]] = i;
	}
	
	void get_row()
	{
		std::string line;
        std::getline(file, line);
        std::vector<std::string> entries(selected_cols.size(), "");
        
        boost::tokenizer<boost::escaped_list_separator<char> > tk(line, boost::escaped_list_separator<char>('\\', colsep, '\"'));
        
        unsigned int n = 0;
        auto tk_it = tk.begin();
        auto scatter_it = scatter_cols.begin();
        for (; tk_it != tk.end(); ++tk_it, ++n, ++scatter_it)
        {
			if (*scatter_it != -1) entries[*scatter_it] = *tk_it;
        }
        
	    if (n == scatter_cols.size()) file_table.push_back(entries);	
	}	
	
	public:
	
	Sequencing_Container(const std::string& input_file, const std::vector<std::string> colselector, char _colsep='\t') : 
		file(input_file),
		colsep(_colsep),
		selected_cols(colselector.size())
	{
		std::cout << "Read headerline\n";
		get_headerline();
		
		auto col_name   = colselector.begin();
		auto col_number = selected_cols.begin();
				
		for (; col_name != colselector.end(); ++col_name, ++col_number)
		{
			auto found_it = std::find_if(headerline.begin(), headerline.end(), [&](const std::string& s1) { return s1.compare(*col_name) == 0; });
			if (found_it == headerline.end()) throw new std::invalid_argument("Cannot find column " + (*col_name) + " within csv file\n");
			*col_number = found_it - headerline.begin();
		}
		
		std::cout << "Selected Cols numbers \n";
		for (auto x=selected_cols.begin(); x != selected_cols.end(); ++x) std::cout << *x << "\t";
		std::cout << std::endl;
		
		invert_scatter();
		
		while (!file.eof()) get_row();
				
		std::cout << "\nFirst 10 interesting rows\n";
		for (auto x=colselector.begin(); x != colselector.end(); ++x) std::cout << *x << "\t";
		std::cout.flush();
		for (auto x=file_table.begin(); x != file_table.end() && x != file_table.begin()+10; ++x) 
		{
			for (auto y = x->begin(); y != x->end(); ++y) std::cout << *y << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	};
	
};
