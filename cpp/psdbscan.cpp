
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <ctime>
#include <ratio>
#include <chrono>
#include <iomanip>
#include <string>
#include <set>
#include <exception>
#include <iostream>

#include "ConfigHolder.hpp"
#include "GeneReader.hpp"
#include "TableHolder.hpp"
#include "BenchmarkSingleton.hpp"


namespace po = boost::program_options;

int main(int argc, char** argv) 
{
	BenchmarkSingleton::Instance().StartLogging();
	BenchmarkSingleton::Instance().StartTimer("Overall");
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("input", po::value<std::vector<std::string> >()->multitoken(), "set input file (csv format, settings may need to be adjusted in config file)")
	("config_file", po::value<std::string>(), "set config file (boost property_tree info file data format)")
	("verbose", "Print some status information during data processing.")
	;
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    
	
	if (vm.count("help")) {
	    std::cout << desc << "\n";
	    return 1;
	}
	
	if (vm.count("verbose")) {
	    BenchmarkSingleton::Instance().BeVerbose();
	}
	
	std::string conf_filename;
	
	if (vm.count("config_file")) 
	{
		conf_filename = vm["config_file"].as<std::string>();
	}
	else
	{
		conf_filename = "data/psdbscan.conf";
	}
	
	std::vector<std::string> input_filenames;
			
	if (vm.count("input")) 
	{
		auto tempvec = vm["input"].as<std::vector<std::string> >();
		input_filenames.assign(tempvec.begin(), tempvec.end());
		
	}
	else
	{
		input_filenames.push_back("data/ExampleData.csv");
	}
	
	
	std::cout << "\n\n***************************************** PSDBSCAN ************************************************\n";
	
	ConfigHolderClass ConfigHolder(conf_filename);
	
	std::for_each(input_filenames.begin(), input_filenames.end(), [&](const std::string& inp_name) { BenchmarkSingleton::Instance().Message("Input file(s) " + inp_name); });
	
	BenchmarkSingleton::Instance().Message("Using configuration file " + conf_filename);
		
	std::cout << "\nReading files ... \n";
	
	BenchmarkSingleton::Instance().StartTimer("Reading Files");
	TableHolderClass TableHolder(input_filenames, ConfigHolder);
	BenchmarkSingleton::Instance().StopTimer("Reading Files");

	std::cout << "\nStarting clustering process... \n";
    BenchmarkSingleton::Instance().StartTimer("Do clustering");
	TableHolder.do_the_clustering();
	BenchmarkSingleton::Instance().StartTimer("Do clustering");
	
	std::cout << "\nWriting output ... \n";
	BenchmarkSingleton::Instance().StartTimer("Writing Output");
	TableHolder.write_to_file();
	BenchmarkSingleton::Instance().StartTimer("Writing Output");
	std::cout << "\nFinished!\n\n";
	BenchmarkSingleton::Instance().StopTimer("Overall");
	std::cout << "***************************************************************************************************\n";
	BenchmarkSingleton::Instance().write_map(ConfigHolder.create_filename(input_filenames[0], ConfigHolderClass::fsBenchmark));
	
	return 0;	
}
