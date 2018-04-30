/* Copyright (c) 2012 Guillaume Collet
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifndef __OPTIONS_H_INCLUDED__
#define __OPTIONS_H_INCLUDED__

#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <stdexcept>

/* This class is a virtual interface for the arguments */
class Arg
{
	protected:
		std::string _small_flag;  // The small version of the argument flag
		std::string _long_flag;   // The long  version of the argument flag
		std::string _description; // The description of the argument
		bool   _needValue;   // This boolean value is true when a value is needed
		bool   _isSet;       // This boolean value is true if the value is setted

		// Protected constructor, only called by inherited classes
		Arg(
			const std::string & small_flag,
			const std::string & long_flag,
			const std::string & description,
			const bool & needed,
			const bool & set
			):
			_small_flag(small_flag),
			_long_flag(long_flag),
			_description(description),
			_needValue(needed),
			_isSet(set)
		{};

	public:

		// The destructor is virtual
		virtual ~Arg(){};

		// Getters
	std::string getSmallFlag() const {
		return _small_flag;
	}
	;
	std::string getLongFlag() const {
		return _long_flag;
	}
	;
	std::string getDescription() const {
		return _description;
	}
	;
		bool   isNeeded()       const {return _needValue;};
		bool   isSetted()       const {return _isSet;};

	/* Try to find the flag in the command_line given in argument,
	 * if the argument is needed but not found in the command_line,
	 * then an runtime_error exception is returned */
	virtual void find(std::vector<std::string> & command_line) = 0;
};


/* ValueArg is an argument with a value of type T */
template <typename T>
class ValueArg : public Arg
{
	private:
		T _value;

	public:
		// The constructor with a default value
		ValueArg(
			const std::string & small_flag,
			const std::string & long_flag,
			const std::string & description,
			const T value
		):
		Arg(small_flag, long_flag, description, false, true),
		_value(value)
		{};

		// The constructor without a default value
		ValueArg(
			const std::string & small_flag,
			const std::string & long_flag,
			const std::string & description
		):
		Arg(small_flag, long_flag, description, true, false)
		{};

		// The destructor does nothing
		virtual ~ValueArg(){};

		// The getter
		T getValue() const {return _value;};

		// The implementation of the "find" virtual function
		void find(std::vector<std::string> & command_line)
		{
			std::vector<std::string>::iterator it = command_line.begin();
			while (it != command_line.end()){
				if (*it == _small_flag || *it == _long_flag){
					it = command_line.erase(it);
					if (it == command_line.end())
					throw std::runtime_error("Value of argument " + _small_flag + ", " + _long_flag + " is missing\n");
				std::istringstream is(*it);
					while (is.good())
						is >> _value;
					command_line.erase(it);
					_isSet = true;
					return;
				}
				it++;
			}
			if (_needValue && !_isSet)
			throw std::runtime_error("Argument " + _small_flag + ", " + _long_flag + " is needed\n");
		}
};


/* SwitchArg is an argument with a boolean (intrinsically not needed)*/
class SwitchArg : public Arg
{
	private:
		bool _value;

	public:
		SwitchArg(
				const std::string & small_flag,
				const std::string & long_flag,
				const std::string & description,
				const bool   & value
			):
			Arg(small_flag,long_flag,description,false,false),
			_value(value)
		{};

		// The getter
		bool getValue() const {return _value;};

		// The implementation of the "find" virtual function
	void find(std::vector<std::string> & command_line)
		{
		std::vector<std::string>::iterator it = command_line.begin();
			while (it != command_line.end()){
				if (*it == _small_flag || *it == _long_flag){
					if (_small_flag == "-h"){
					throw std::runtime_error("");
					}
					command_line.erase(it);
					_value = true;
					return;
				}
				it++;
			}
		}
};

/*
 * The Options class manages the argument given on the command line.
 * The implementation is static and then accessible from anywhere.
 */
class Options
{
	private:
	std::string appName;
	std::vector<std::string> command_line;
	std::map<std::string, Arg*> arg_list; // Map of argument objects sorted by their flag

		// The constructor is private
		Options():
			input_fname  (""),
			matrix_fname (""),
			output_fname ("output.txt"),
			statistic    ("wentropy"),
			nb_seq       (500),
			verbose      (false),
			global       (false),
			threshold    (0.8),
			factor_a     (1.0),
			factor_b     (0.5),
			factor_c     (3.0),
			window       (3)
		{};

		// The accessor to the uniq Options object
		static Options & GetNC()
		{
			static Options opt;
			return opt;
		}

		// getEnvVar : Get an environment variable of name env
	std::string getEnvVar(std::string env)
		{
			char * env_p;
		std::string env_s;
			env_p = getenv(env.c_str());
			if (env_p != NULL)
				env_s = env_p;
			else
			std::cerr << "Warning: Environment variable " << env << " is not found\n";
			return env_s;
		}

		// Reduce a pathname in a basename
	std::string basename(std::string fname)
		{
			int pos = (int) fname.find_last_of('/');
			return fname.substr(pos+1, fname.size() - pos);
		};

		// Initialize the options (This is where we add command line arguments)
		void Init(int argc, char * const argv[])
		{
			try {
				// Set the command_line object
				for (int i = 1; i < argc; i++){
				command_line.push_back(std::string(argv[i]));
				}
				// Set the application name
				appName = basename(argv[0]);

				// Get the environment variable to find the scoring matrix
			std::string smat_path = getEnvVar("SCORE_MAT_PATH");
				if (smat_path.empty()) {
					smat_path = "data/aaindex";
				std::cerr << "SCORE_MAT_PATH is not set, use data/aaindex instead\n";
				}
				/*
				 * 2 sorts of arguments can be added:
				 *   - A ValueArg  which is for flags with values
				 *   - A SwitchArg which is for boolean flags
				 *
				 * Be careful with ValueArg declaration, the fourth argument is optional
				 * if given, it will be the default value for the argument.
				 */

				//1 - create the argument as a ValueArg or SwitchArg.
			ValueArg<std::string> iArg("-i", "--input", "MSA input file name");
			ValueArg<std::string> mArg("-m", "--matrix", "Score matrix file name", smat_path + "/HENS920102.mat");
			ValueArg<std::string> oArg("-o", "--output", "Output file name [default=ouput.txt]", "output.txt");
			ValueArg<std::string> sArg("-s", "--statistic", "Statistics [default=wentropy]", "wentropy");
				ValueArg<int>    nArg("-n", "--nb_seq",    "Maximum number of sequences read [default=500]",     500);
				SwitchArg        vArg("-v", "--verbose",   "Verbose mode",                                     false);
				SwitchArg        gArg("-g", "--global",    "Output the global score",                          false);
				SwitchArg        hArg("-h", "--help",      "Print this help",                                  false);
				ValueArg<float>  tArg("-t", "--threshold", "Threshold to print correlation [default=0.8]",       0.8);
				ValueArg<float>  aArg("-a", "--trident_a", "Factor applied to t(x) (see trident) [default=1.0]", 1.0);
				ValueArg<float>  bArg("-b", "--trident_b", "Factor applied to r(x) (see trident) [default=0.5]", 0.5);
				ValueArg<float>  cArg("-c", "--trident_c", "Factor applied to g(x) (see trident) [default=3.0]", 3.0);
				ValueArg<int>    wArg("-w", "--window",    "Number of side columns (jensen score)",                3);

				// 2 -  add the argument to the arg_list for further use (print_usage).
				arg_list[iArg.getSmallFlag()] = &iArg;
				arg_list[mArg.getSmallFlag()] = &mArg;
				arg_list[oArg.getSmallFlag()] = &oArg;
				arg_list[sArg.getSmallFlag()] = &sArg;
				arg_list[nArg.getSmallFlag()] = &nArg;
				arg_list[vArg.getSmallFlag()] = &vArg;
				arg_list[gArg.getSmallFlag()] = &gArg;
				arg_list[hArg.getSmallFlag()] = &hArg;
				arg_list[tArg.getSmallFlag()] = &tArg;
				arg_list[aArg.getSmallFlag()] = &aArg;
				arg_list[bArg.getSmallFlag()] = &bArg;
				arg_list[cArg.getSmallFlag()] = &cArg;
				arg_list[wArg.getSmallFlag()] = &wArg;

				// 3 - try to find the argument in the command line to set up the value.
				hArg.find(command_line);
				iArg.find(command_line);
				mArg.find(command_line);
				oArg.find(command_line);
				sArg.find(command_line);
				nArg.find(command_line);
				vArg.find(command_line);
				gArg.find(command_line);
				tArg.find(command_line);
				aArg.find(command_line);
				bArg.find(command_line);
				cArg.find(command_line);
				wArg.find(command_line);

				// If something is left in the command line... It is not an argument of the program -> error
				if (command_line.size() > 0){
				throw std::runtime_error("Unknown flag: " + command_line[0] + "\n");
				}

				// Get arguments
				input_fname  = iArg.getValue();
				matrix_fname = mArg.getValue();
				output_fname = oArg.getValue();
				statistic    = sArg.getValue();
				nb_seq       = nArg.getValue();
				verbose      = vArg.getValue();
				global       = gArg.getValue();
				threshold    = tArg.getValue();
				factor_a     = aArg.getValue();
				factor_b     = bArg.getValue();
				factor_c     = cArg.getValue();
				window       = wArg.getValue();
		} catch (std::exception &e) {
				throw;
			}
		}

	public:
		/* List of options */
	std::string input_fname;  // The file name of the multiple alignment */
	std::string matrix_fname; // The file name of the scoring matrix */
	std::string output_fname; // The name of the output file */
	std::string statistic;    // The name of the statistic */
		int    nb_seq;       // The number of sequences to read in the multiple alignment */
		bool   verbose;      // The switch for verbose mode */
		bool   global;       // The switch to output only the global alignment score */
		float  threshold;    // The threshold for correlation print */
		float  factor_a;     // The factor applied to the first  member of trident score */
		float  factor_b;     // The factor applied to the second member of trident score */
		float  factor_c;     // The factor applied to the third  member of trident score */
		int    window;       // The size of the window to take in account side columns (jensen stat only) */

		/* Universal accessor */
		static Options const & Get()
		{
			return GetNC();
		}

		/* Initialization */
		static void Parse(int argc_, char * const argv_[])
		{
			Options & opt = GetNC();
			try{
				opt.Init(argc_,argv_);
		} catch (std::exception &e) {
				throw;
			}
		}

		/* Print usage */
		static void print_usage()
		{
			Options & opt = GetNC();
		std::map<std::string, Arg*>::iterator it = opt.arg_list.begin();
			int sflag_size = 0;
			int lflag_size = 0;
			int desc_size  = 0;
			while (it != opt.arg_list.end()){
				if (sflag_size < (int) it->second->getSmallFlag().length()) {
					sflag_size = (int) it->second->getSmallFlag().size();
				}
				if (lflag_size < (int) it->second->getLongFlag().size() ) {
					lflag_size = (int) it->second->getLongFlag().size();
				}
				if (desc_size < (int) it->second->getDescription().size()) {
					desc_size = (int) it->second->getDescription().size();
				}
				it++;
			}
		std::cerr << "Usage: " << opt.appName;
			it = opt.arg_list.begin();
			int nb = 1;
			while (it != opt.arg_list.end()){
				if (it->second->isNeeded()) {
				std::cerr << " " << it->second->getSmallFlag() << " file";
					nb++;
				}
				it++;
			}
		std::cerr << " [options]\n\n";
		std::cerr << "Available statistics: \n";
		std::cerr << "  sumofpairs (1)\n";
		std::cerr << "  wentropy (1)\n";
		std::cerr << "  trident  (1)\n";
		std::cerr << "  mvector  (1)\n";
		std::cerr << "  jensen   (1)\n";
		std::cerr << "  kabat    (1)\n";
		std::cerr << "  gap      (1)\n";
		std::cerr << "\nOptions:\n";
			it = opt.arg_list.begin();
			while (it != opt.arg_list.end()){
			std::string flag = it->second->getSmallFlag() + ",";
			std::cerr << "   " << std::setw(sflag_size + 1) << std::left << flag;
			std::cerr << " " << std::setw(lflag_size) << std::left << it->second->getLongFlag();
			std::cerr << " : " << std::setw(desc_size) << std::left << it->second->getDescription();
			std::cerr << "\n";
				it++;
			}
		std::cerr << "\n";
		}
};

#endif
