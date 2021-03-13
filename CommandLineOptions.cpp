#include "CommandLineOptions.h"
#include <iostream>
#include <cstdlib>

using namespace std;

#include <boost/program_options.hpp>

namespace po = boost::program_options;

bool Options(int argc, char** argv, po::variables_map &vm, string &InitialCondition, double &dt, double &T, double &h){
	
	try{
		
		po::options_description desc("Available Options");
		desc.add_options()
			("help", "Print help message")
			("ic-dam-break", "Use dam-break initial condition")
			("ic-droplet", "Use droplet initial condition")
			("ic-block-drop", "Use block drop initial condition")
			("ic-one-particle", "Use one particle validation case")
			("ic-two-particles", "Use two particle validation case")
			("ic-four-particles", "Use four particle validation case")
			("dt", po::value<double>(), "Select time step")
			("T", po::value<double>(), "Total integration time")
			("h", po::value<double>(), "Radius of influence per particle");
	
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	
		if (vm.count("help")){
			cout << desc << endl;
			return true;
		}
		
		if (vm.count("ic-dam-break")){
			InitialCondition = "ic-dam-break";
		}
		
		if (vm.count("ic-droplet")){
			InitialCondition = "ic-droplet";
		}
		
		if (vm.count("ic-block-drop")){
			InitialCondition = "ic-block-drop";
		}
		
		if (vm.count("ic-one-particle")){
			InitialCondition = "ic-one-particle";
		}
		
		if (vm.count("ic-two-particles")){
			InitialCondition = "ic-two-particles";
		}
		if (vm.count("ic-four-particles")){
			InitialCondition = "ic-four-particles";
		} 
		
		if ((vm["dt"].as<double>() <= 0) || (vm["T"].as<double>() <= 0) || (vm["h"].as<double>() <= 0)){
			cout << "Please provide positive values of dt, T and h. " << endl;
			return false;
		}
		else {
			dt = vm["dt"].as<double>();
			T  = vm["T"].as<double>();
			h  = vm["h"].as<double>();
		}
		
	}
	//Throw exceptions if there are invalid inputs
	catch(exception const &e){
		
		cout << e.what() << endl;
		return false;
	}
	
	return true;
	
	
}



