/* Copyright (c) 2010 Guillaume Collet
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

using namespace std;

class Arg {
protected:
	string	_small_flag;	/**< The small version of the argument flag */
	string	_long_flag;		/**< The long version of the argument flag */
	string	_description;	/**< The description of the argument (for usage printing)*/
	bool		_needValue;		/**< This boolean value is true when a value is needed */
	bool		_isSet;
	/* Protected constructor, only called by inherited classes */
	Arg(const string & small_flag, const string & long_flag, const string & description, const bool & needed, const bool & set):
	_small_flag(small_flag), _long_flag(long_flag),_description(description),	_needValue(needed), _isSet(set)
	{};
	
public:
	/**
	 * The constructor does nothing
	 */
	Arg(){};
	
	
	/**
	 * The destructor is virtual
	 */
	virtual ~Arg(){};
	
	/**
	 * Returns the argument small flag.
	 */
	const string getSmallFlag() const {return _small_flag;};
	
	/**
	 * Returns the argument long flag.
	 */
	const string getLongFlag() const {return _long_flag;};
	
	/**
	 * Returns the argument description.
	 */
	string getDescription() const {return _description;};
	
	/**
	 * Indicates whether the argument needs a value or not.
	 */
	bool isNeeded() const {return _needValue;};
	
	/**
	 * Indicates whether the argument is setted or not.
	 */
	bool isSetted() const {return _isSet;};
	
	/**
	 * Sets the value of the argument.
	 */
	virtual void setValue(string val){};
	
	/**
	 * Try to find the flag in the command line
	 * if needed and not found -> error
	 */
	virtual void find(vector<string> & command_line){};
};

/**
 * ValueArg is an argument with a value of type T
 */
template <typename T>
class ValueArg : public Arg {
private:
	T _value;
	
public:
	ValueArg(const string & small_flag,
					 const string & long_flag,
					 const string & description,
					 const T value
					 ):Arg(small_flag,long_flag,description,false,true), _value(value){};
	
	ValueArg(const string & small_flag,
					 const string & long_flag,
					 const string & description
					 ):Arg(small_flag,long_flag,description,true,false){};
	
	T getValue() const {return _value;};
	
	void setValue(string val) {
		istringstream is(val);
		while (is.good()) {
			is >> _value;
		}
	}
	
	void  find(vector<string> & command_line){
		vector<string>::iterator it = command_line.begin();
		while (it != command_line.end()){
			if (*it == _small_flag || *it == _long_flag){
				it = command_line.erase(it);
				if (it == command_line.end()) {
					throw runtime_error("Value of argument " + _small_flag + ", " + _long_flag + " is missing\n");
				}
				istringstream is(*it);
				while (is.good()) {
					is >> _value;
				}
				command_line.erase(it);
				_isSet = true;
				return;
			}
			it++;
		}
		if (_needValue && !_isSet){
			throw runtime_error("Argument " + _small_flag + ", " + _long_flag + " is needed\n");
		}
	}
};


/**
 * SwitchArg is an argument with a boolean 
 * it is intrinsically not needed -> _needValue = false
 */
class SwitchArg : public Arg {
private:
	bool _value;
public:
	SwitchArg(const string & small_flag,
						const string & long_flag,
						const string & description,
						const bool & value):
	Arg(small_flag,long_flag,description,false,false),
	_value(value)
	{};
	
	bool getValue() const {return _value;};
	void setValue(string val) {_value = true;};
	
	void  find(vector<string> & command_line){
		vector<string>::iterator it = command_line.begin();
		while (it != command_line.end()){
			if (*it == _small_flag || *it == _long_flag){
				if (_small_flag == "-h"){
					throw runtime_error("");
				}
				command_line.erase(it);
				_value = true;
				return;
			}
			it++;
		}
	}
	
};

class Options
{
private:
	string appName;
	vector<string> command_line;
	map<string, Arg> arg_list;	/**< Map of argument objects by their flag */
	
	
	/*
	 * Constructor is private
	 */
	Options(){};
  
	/* 
	 * Accessor for the uniq Options object
	 */
	static Options & GetNC(){
		static Options opt;
		return opt;
	}
  
	/*
	 * getEnvVar : Get an environment variable of name env
	 */
	string getEnvVar(string env){
		char * env_p;
		string env_s;
		env_p = getenv(env.c_str());
		if (env_p != NULL){
			env_s = env_p;
		} else {
			throw runtime_error("Error: Environment variable " + env + " is needed\n");
		}
		return env_s;
	}
    
	/*
	 * Reduce a pathname in a basename
	 */
	string basename(string fname){
		int pos = fname.find_last_of('/');
		return fname.substr(pos+1, fname.size() - pos);
	};
	
	/*
   * Initialize the options (This is where you edit the file)
   */
	void Init(int argc, char * const argv[]) {
		try {
			
			// Set the command_line object
			for (int i = 1; i < argc; i++){
				command_line.push_back(string(argv[i]));
			}
			// Set the application name
			appName = basename(argv[0]);
			
			/*
			 * You can add 2 kind of arguments :
			 *   A ValueArg which is for flags with values
			 *   A SwitchArg which is for boolean flags
			 *
			 * Be careful with ValueArg declaration, the 3 first parameters must be there: small_flag, long_flag, description.
			 * The fourth is optional and is a default value for the argument. 
			 * If it is not given, then the flag must be in the command line to set the value (-i option in this example).
			 */
			
			/** Set options
			 *	3 lines are needed to set up an option
			 */
			
			//1 - create the argument as a ValueArg or SwitchArg.
			ValueArg<string>	smArg("-sm",	"--score_matrix", "Score matrix file", "mclachlan71.mat");
			ValueArg<string>	smpArg("-sp",	"--score_matrix_path", "Path to the score matrix", getEnvVar("SCORE_MAT_PATH"));
			ValueArg<string>	mafArg("-m",	"--mult_ali", "File name of the multiple alignment");
			ValueArg<string>	outArg("-o",	"--output", "Name of the output file [default=mali_file]", "output.txt");
			ValueArg<string>	msArg("-s",		"--statistic", "Name of the statistic used [default=trident]", "wentropy");
			ValueArg<float>		thArg( "-t",	"--threshold", "Threshold to print correlation [default=0.8]", 0.8);
			ValueArg<float>		taArg( "-a",	"--trident_a", "Factor applied to t(x) (see trident) [default=1.0]", 1.0);
			ValueArg<float>		tbArg( "-b",	"--trident_b", "Factor applied to r(x) (see trident) [default=0.5]", 0.5);
			ValueArg<float>		tcArg( "-c",	"--trident_c", "Factor applied to g(x) (see trident) [default=3.0]", 3.0);
			SwitchArg					vArg(  "-v",	"--verbose", "Verbose mode", false);
			SwitchArg					gArg(  "-g",	"--global", "Output only the global score of the multiple alignment", false);
			SwitchArg					bArg(  "-ba",	"--basic", "Add a basic output nbAA/col", false);
			SwitchArg					hArg(  "-h",	"--help", "Print this help", false);
			ValueArg<int>			nsArg( "-n",	"--nb_seq", "Maximum number of sequences read [default=500]", 500);
			ValueArg<int>			wArg(  "-w",	"--window", "Number of side columns (jensen score)", 3);
			
			// 2 -  add the argument to the arg_list for further use (print_usage).
			arg_list[smArg.getSmallFlag()]  = smArg;
			arg_list[smpArg.getSmallFlag()] = smpArg;
			arg_list[mafArg.getSmallFlag()] = mafArg;
			arg_list[outArg.getSmallFlag()] = outArg;
			arg_list[msArg.getSmallFlag()]  = msArg;
			arg_list[thArg.getSmallFlag()]  = thArg;
			arg_list[taArg.getSmallFlag()]  = taArg;
			arg_list[tbArg.getSmallFlag()]  = tbArg;
			arg_list[tcArg.getSmallFlag()]  = tcArg;
			arg_list[vArg.getSmallFlag()]   = vArg;
			arg_list[gArg.getSmallFlag()]   = gArg;
			arg_list[bArg.getSmallFlag()]   = bArg;
			arg_list[nsArg.getSmallFlag()]  = nsArg;
			arg_list[wArg.getSmallFlag()]   = wArg;
			arg_list[hArg.getSmallFlag()]   = hArg;
			
			// 3 - try to find the argument in the command line to set up the value.
			hArg.find(command_line);
			mafArg.find(command_line);
			smArg.find(command_line);
			smpArg.find(command_line);
			outArg.find(command_line);
			msArg.find(command_line);
			thArg.find(command_line);
			taArg.find(command_line);
			tbArg.find(command_line);
			tcArg.find(command_line);
			vArg.find(command_line);
			gArg.find(command_line);
			bArg.find(command_line);
			nsArg.find(command_line);
			wArg.find(command_line);
			
			// If something is left in the command line... It is not an argument of the program -> error
			if (command_line.size() > 0){
				throw runtime_error("Unknown flag: " + command_line[0] + "\n");
			}
			
			// Get arguments
			score_matrix_fname  = smArg.getValue();
			score_matrix_path   = smpArg.getValue();
			mult_ali_fname			= mafArg.getValue();
			output_name         = outArg.getValue();
			statistic						= msArg.getValue();
			threshold		        = thArg.getValue();
			factor_a						= taArg.getValue();
			factor_b						= tbArg.getValue();
			factor_c						= tcArg.getValue();
			verbose							= vArg.getValue();
			global							= gArg.getValue();
			basic								= bArg.getValue();
			nb_seq							= nsArg.getValue();
			window							=	wArg.getValue();
			
			/* Changes the output name if not given by the user */
			if (output_name == "output.txt"){
				output_name = mult_ali_fname.substr(0,mult_ali_fname.find("."));
				output_name += ".stat";
			}
			
			
		} catch (exception &e) {
        throw;
		}
	}

public:
	/*
	 * List of options (Edit to add options)
	 */
	/* List of options */
	string score_matrix_fname;  /**< The file name of the scoring matrix */
	string score_matrix_path;  	/**< The path to the scoring matrix */
	string mult_ali_fname;			/**< The file name of the multiple alignment */
	string output_name;         /**< The name of the output file */
	string statistic;           /**< The name of the statistic */
	float threshold;            /**< The threshold for correlation print */
	float factor_a;							/**< The factor applied to the first  member of trident score */
	float factor_b;             /**< The factor applied to the second member of trident score */
	float factor_c;							/**< The factor applied to the third  member of trident score */
	bool verbose;								/**< The switch for verbose mode */
	bool global;								/**< The switch to output only the global alignment score */
	bool basic;									/**< The switch for basic information output */
	int nb_seq;									/**< The number of sequences to read in the multiple alignment */
	int window;                 /**< The size of the window to take in account side columns (jensen stat only) */
	
	/* Universal accessor */
	static Options const & Get(){
		return GetNC();
	}
  
	/* Initialization */
	static void Parse(int argc_, char * const argv_[]) {
		Options & opt = GetNC();
		try{
			opt.Init(argc_,argv_);
		} catch (exception &e) {
			throw;
		}
	}
	
	/* Print usage */
	static void print_usage(){
		Options & opt = GetNC();
		
		map<string,Arg>::iterator it = opt.arg_list.begin();
		int sflag_size = 0;
		int lflag_size = 0;
		int desc_size = 0;
		while (it != opt.arg_list.end()){
			if (it->second.getSmallFlag().size() > sflag_size) {
				sflag_size = it->second.getSmallFlag().size();
			}
			if (it->second.getLongFlag().size() > lflag_size) {
				lflag_size = it->second.getLongFlag().size();
			}
			if (it->second.getDescription().size() > desc_size) {
				desc_size = it->second.getDescription().size();
			}
			it++;
		}
		cerr << "Usage: " << opt.appName ;
		it = opt.arg_list.begin();
		int nb = 1;
		while (it != opt.arg_list.end()){
			if (it->second.isNeeded()) {
			  cerr << " " << it->second.getSmallFlag() << " file";
				nb++;
			}
			it++;
		}
		
		cerr << " [options]\n\n";
		cerr << "Available statistics: \n";
		cerr << "  wentropy (1)\n";
		cerr << "  trident  (1)\n";
		cerr << "  mvector  (1)\n";
		cerr << "  jensen   (1)\n";
		cerr << "  kabat    (1)\n";
		cerr << "  gap      (1)\n";
		cerr << "\nOptions:\n";
		it = opt.arg_list.begin();
		while (it != opt.arg_list.end()){
			string flag = it->second.getSmallFlag() + ",";
			cerr << "   " << setw(sflag_size+1) << left << flag;
			cerr << " " << setw(lflag_size) << left << it->second.getLongFlag()<< " : ";
			cerr << setw(desc_size) << left << it->second.getDescription() << "\n";
			it++;
		}
		cerr << "\n";
	}
};

#endif

