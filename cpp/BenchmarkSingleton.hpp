/*
 * BenchmarkSingleton.hpp
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
 #include <fstream>
 #include <chrono>
 #include <map>
 #include <string>
 #include <tuple>
 
class BenchmarkSingleton
{
	using timer_value_type = std::pair<std::chrono::high_resolution_clock::time_point, std::chrono::high_resolution_clock::time_point>;
		
	std::map <std::string, timer_value_type> timermap;
	std::vector<std::pair<std::string, std::string> > data;
		
	bool	verbose;
		
	BenchmarkSingleton() : verbose(false) {}
	
	public:
	static BenchmarkSingleton& Instance()
	{
		static BenchmarkSingleton obj;
		return obj;
	};
	
	void StartTimer(std::string&& TimerName)
	{
		auto pos = timermap.find(TimerName);
		if (pos != timermap.end())
		{
			(*pos).second.first = std::chrono::high_resolution_clock::now();
		} 
		else
		{
			timermap.emplace(TimerName, std::make_pair(std::chrono::high_resolution_clock::now(),std::chrono::high_resolution_clock::now()));
		}
	}
	
	void StopTimer(std::string&& TimerName)
	{
		auto pos = timermap.find(TimerName);
		if (pos != timermap.end())
		{
			(*pos).second.second = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> tdiff = std::chrono::duration_cast<std::chrono::duration<double>>((*pos).second.second - (*pos).second.first);
			
			data.emplace_back(std::make_pair(TimerName, std::to_string(tdiff.count()) + " sec"));			
			if (verbose) std::cout << "\nTimer " << TimerName << ": " << tdiff.count() << " sec";
		}
	}
	
	void StartLogging()
	{
		StartTimer("Logging");
	}
	
	void BeVerbose()
	{
		verbose = true;
	}
	
	void Message(std::string&& Message)
	{
		auto pos = timermap.find("Logging");
		
		if (pos == timermap.end()) {
			StartTimer("Logging");
			pos = timermap.find("Logging");
		}
		
		(*pos).second.second = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> tdiff = std::chrono::duration_cast<std::chrono::duration<double>>((*pos).second.second - (*pos).second.first);
		
		data.emplace_back( std::make_pair("At " + std::to_string(tdiff.count()) + " sec", Message));
		if (verbose) std::cout << "\n" << Message;
	}
	
	
	void write_map(const std::string& filename)
	{
		std::ofstream out(filename, std::ofstream::out);
		// header
		
		out << "Timer/Message" << "\t" << "Value" << "\n";
		for (auto it = data.begin(); it != data.end(); ++it)
		{
			out << (*it).first << "\t" << (*it).second << "\n";
		}
		out.close();
	}
	
};
