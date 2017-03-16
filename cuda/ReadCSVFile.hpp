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
#include <queue>
#include <algorithm>
#include <thread>
#include <condition_variable>
#include <mutex>


class Sequencing_Container
{
	std::vector<std::vector <std::string>> file_table;
	std::vector<std::string> headerline, final_col_names;
	
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
	
	
	typedef std::pair<std::string, unsigned int> waiting_queue_element_t;
	std::queue<waiting_queue_element_t> waiting_queue;
	
	std::condition_variable cv, cv2;
	std::mutex queue_mutex, write_mutex;
	bool finished_reading;
	unsigned current_entry;
	
	
	void push_unprocessed_row(waiting_queue_element_t entry)
	{
		std::unique_lock<std::mutex> lck(queue_mutex);
		waiting_queue.push(entry);
	    cv.notify_all();
	}
	
	void process_entries()
	{
		while (!finished_reading || !waiting_queue.empty())  // tasks to do?
		{
			// wait for lock and task
			std::unique_lock<std::mutex> lck(queue_mutex);
			
	        while (!finished_reading && waiting_queue.empty()) cv.wait(lck);
	        
	        // take task and release lock
	        if (!waiting_queue.empty())
	        {
		        waiting_queue_element_t xentry = waiting_queue.front();
		        waiting_queue.pop();
		        cv.notify_all();
		        lck.unlock();
		        
		        // process task
		        std::vector<std::string> entries(selected_cols.size(), "");
	            boost::tokenizer<boost::escaped_list_separator<char> > tk(xentry.first, boost::escaped_list_separator<char>('\\', colsep, '\"'));
	            unsigned int n = 0;
				auto tk_it = tk.begin();
				auto scatter_it = scatter_cols.begin();
				for (; tk_it != tk.end(); ++tk_it, ++n, ++scatter_it)
				{
					if (*scatter_it != -1) entries[*scatter_it] = *tk_it;
				}
	        
	            // try to write line
				std::unique_lock<std::mutex> lck2(write_mutex);
				while (xentry.second > current_entry) cv2.wait(lck2);
				
				if (n == scatter_cols.size()) 
				{
					auto ft_it    = file_table.begin();
					auto entry_it = entries.begin();
					for (; ft_it != file_table.end(); ++ft_it, ++entry_it) ft_it->push_back(*entry_it);
				} 
				++current_entry;
		        cv2.notify_all();
			}
		}
	}
	
	
	
	void produce_rows()
	{
		std::string line;
		finished_reading = false;  // 1st initialization
		current_entry = 0; // 1st initialization
		unsigned int line_number = 0;
		
		// start threads
		std::vector<std::thread> xthreads;
		for (int i=0; i<std::thread::hardware_concurrency(); ++i) xthreads.emplace_back(&Sequencing_Container::process_entries,this);

		while (!file.eof()) {
			std::getline(file, line);
			push_unprocessed_row(std::make_pair(line, line_number));
			++line_number;
		}
		
		std::unique_lock<std::mutex> lck(queue_mutex);
		finished_reading = true;
		cv.notify_all();
		lck.unlock();
		
		// wait for threads to finish
		
		for (auto& x: xthreads) x.join();
	}
    
	public:
	
	Sequencing_Container(const std::string& input_file, const std::vector<std::string> colselector, char _colsep='\t') : 
		file(input_file),
		colsep(_colsep),
		selected_cols(colselector.size()),
		final_col_names(colselector),
		file_table(colselector.size())
	{
		get_headerline();
		
		auto col_name   = colselector.begin();
		auto col_number = selected_cols.begin();
				
		for (; col_name != colselector.end(); ++col_name, ++col_number)
		{
			auto found_it = std::find_if(headerline.begin(), headerline.end(), [&](const std::string& s1) { return s1.compare(*col_name) == 0; });
			if (found_it == headerline.end()) throw new std::invalid_argument("Cannot find column " + (*col_name) + " within csv file\n");
			*col_number = found_it - headerline.begin();
		}
		
		
		invert_scatter();
		
		produce_rows();
		
		/*		
		std::cout << "\nFirst 10 interesting rows\n";
		for (auto x=final_col_names.begin(); x != final_col_names.end(); ++x) std::cout << std::left << std::setw(30) << *x << "\t";
		std::cout << std::endl; std::cout.flush();
		for (unsigned i = 0; i < 10; ++i)
		{
			for (auto ft_it = file_table.begin(); ft_it != file_table.end(); ++ft_it)
				std::cout << std::left << std::setw(30) << (*ft_it)[i] << "\t";
			std::cout << std::endl;
		}
		std::cout << "******";
		std::cout.width(0);
		std::cout << std::endl;
		*/
	};
	
	std::vector<std::string> get_column(std::string&& which_col)
	{
		int col_index = std::find_if(final_col_names.begin(), final_col_names.end(), [&](std::string& sc) { return sc.compare(which_col) == 0;} ) - final_col_names.begin();
		if (col_index >= final_col_names.size()) throw std::runtime_error("Column name not found");
		return file_table[col_index];
	}
	
};
