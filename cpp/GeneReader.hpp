/*
 * GeneReader.hpp
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
#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<stdexcept>

/*
 * The following class reads a simple text file containing gene names
 * for identification of genes in a specified column of an immune 
 * repertoire analysis table like an IMGT output file.
 * 
 * The gene text file contains all gene names with gene names in lines.
 * So the format is
 * 
 * IGHV1-23*
 * IGHV1-24*
 * ...
 * 
 * Comment: the last character in a row should be a star or a blank to
 * have a delimiter! Otherwise, the algorithm cannot distinguish e.g.
 * IGHV1-2 from IGHV1-23!
*/

using namespace std;

class Gene_Reader 
{
	vector<string> gene_names_list;
	ifstream       file;
	bool empty;

	public:

	Gene_Reader() : empty(true)
	{
		// do nothing
	}
	
	

	Gene_Reader(const string& filename) : file(filename), empty(false)
	{
		if(!file.is_open())
		{
		  throw runtime_error("Cannot read from file " + filename);
		}
		while (!file.eof()) {
			string line;
			getline(file, line);
			if (line.compare("") != 0) gene_names_list.push_back(line);
		}
	}
	
	auto begin() -> decltype(gene_names_list.begin())
	{
		return gene_names_list.begin();
	}
	
	auto end() -> decltype(gene_names_list.end())
	{
		return gene_names_list.end();
	}
	

	bool is_empty() 
	{
		return empty;
	}

};

