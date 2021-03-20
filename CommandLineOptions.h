#ifndef COMMANDLINEOPTIONS_H
#define COMMANDLINEOPTIONS_H

using namespace std;

#include <boost/program_options.hpp>

namespace po = boost::program_options;

void Options(int argc, char** argv, po::variables_map &vm, string&InitialCondition, double &dt, double &T, double &h);



#endif // COMMANDLINEOPTIONS_H
