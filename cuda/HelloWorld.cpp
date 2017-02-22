/*
 * HelloWorld.cpp
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
 
#include <iostream>
#include <string>
#include <ctime>
#include <ratio>
#include <chrono>
#include "boost/program_options.hpp"
#include "ReadCSVFile.hpp"

int main(int argc, char** argv) {
	   // Declare the supported options.
	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "produce help message")
	    ("input", boost::program_options::value<std::string>(), "set input file (csv file containing sequencing data)")
	    ("output", boost::program_options::value<std::string>(), "set output file (csv file containing columns of input file together with clone id column)")
	;
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);    
	
	if (vm.count("help")) {
	    std::cout << desc << "\n";
	    return 1;
	}
	
	if (vm.count("input")) {
	    std::cout << "Input file: " << vm["input"].as<std::string>() << "\n";
	    
	    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	    
	    Sequencing_Container SeqCont(vm["input"].as<std::string>(), {"Sequence.number", "Sequence.ID", "V-GENE.and.allele", "J-GENE.and.allele", "D-GENE.and.allele", "CDR3-IMGT"} );
	    
	    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	    
	    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

        std::cout << "It took me " << time_span.count() << " seconds.";
        std::cout << std::endl;
	    std::cout.flush();    
	} else {
	    std::cout << "Input file missing.\n";
	}
	
	if (vm.count("output")) {
	    std::cout << "Output file: " << vm["output"].as<std::string>() << "\n";
	} else {
	    std::cout << "Output file missing, using standard.\n";
	}
   
} 
