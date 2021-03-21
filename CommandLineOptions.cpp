#include "CommandLineOptions.h"
#include <iostream>
#include <cstdlib>

using namespace std;

#include <boost/program_options.hpp>

namespace po = boost::program_options;

void Options(int argc, char** argv, po::variables_map &vm, string &InitialCondition, double &dt, double &T, double &h){
		
		//Initialize dummy count variable to help catch exceptions
		int dummycount = 0;
		
		po::options_description desc("Available Options");
		desc.add_options()
			("help", "Print help message")
			("ic-dam-break", "Use dam-break initial condition")
			("ic-droplet", "Use droplet initial condition")
			("ic-block-drop", "Use block drop initial condition")
			("ic-one-particle", "Use one particle validation case")
			("ic-two-particles", "Use two particles validation case")
			("ic-three-particles", "Use three particles validation case")
			("ic-four-particles", "Use four particles validation case")
			("dt", po::value<double>()->default_value(0.0001), "Select time step")
			("T", po::value<double>()->default_value(3), "Total integration time")
			("h", po::value<double>()->default_value(0.01), "Radius of influence per particle");
	
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	
		if (vm.count("help")){
			
			cout << desc << endl;
			
		}
		
		if (vm.count("ic-dam-break")){
			
			InitialCondition = "ic-dam-break";
			
			dummycount++;
		
		}
		
		if (vm.count("ic-droplet")){
			
			InitialCondition = "ic-droplet";
			
			dummycount++;
		
		}
		
		if (vm.count("ic-block-drop")){
			
			InitialCondition = "ic-block-drop";
	
			dummycount++;
			
		}
		
		if (vm.count("ic-one-particle")){
		
			InitialCondition = "ic-one-particle";
			
			dummycount++;
		
		}
		
		if (vm.count("ic-two-particles")){
			
			InitialCondition = "ic-two-particles";
			
			dummycount++;
		
		}
		if (vm.count("ic-four-particles")){
			
			InitialCondition = "ic-four-particles";
			
			dummycount++;
		
		} 
		
		if ((vm["dt"].as<double>() <= 0) || (vm["T"].as<double>() <= 0) || (vm["h"].as<double>() <= 0)){
			
			throw std::logic_error("dt, T and h must be greater than 0");
			
			
		
		}
		else {
			
			dt = vm["dt"].as<double>();
			
			T  = vm["T"].as<double>();
			
			h  = vm["h"].as<double>();
		
		}
		
		if (dummycount > 1){
			
			throw std::logic_error("Can only have 1 initialization condition.");
			
		}
	
	
}



