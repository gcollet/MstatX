/*
 *  options.h
 *  MstatX++
 *
 *  Created by gcollet on 22/06/10.
 *  Copyright 2010 IRISA. All rights reserved.
 *
 */

#ifndef __OPTIONS_H_INCLUDED__
#define __OPTIONS_H_INCLUDED__

#include <string>
#include <cstdlib>
#include <stdexcept>
#include "arg_line.h"

using namespace std;

class Options
{
  private:
  
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
      }
      return env_s;
    }
    
    
   /*
    * Initialize the options
    */
    void Init(int argc, char * const argv[]) {
      try {
        // Declare the command line
        ARGU::ArgLine argline(argc, argv);
				
        // Set options in ArgLine object
				ARGU::ValueArg<string>	mafArg(	"-ma"			,"--mult_ali"					,"File name of the multiple alignment"								,"string"	,"");
				argline.addArg(&mafArg);
				ARGU::ValueArg<string>	smpArg(	"-sp"			,"--score_matrix_path","Path to the score matrix"														,"string" ,getEnvVar("SCORE_MAT_PATH"));
				argline.addArg(&smpArg);
				ARGU::SwitchArg					verbArg("-v"			,"--verbose"					,"Verbose mode"																				,"string");
				argline.addArg(&verbArg);
				ARGU::ValueArg<string>	outArg(	"-o"			,"--output"						,"Name of the output file [default=mali_file]"				,"string"	,"");
				argline.addArg(&outArg);
				ARGU::ValueArg<string>	msArg(	"-s"			,"--statistic"				,"Name of the statistic used [default=mlc]"						,"string"	,"trident");
				argline.addArg(&msArg);
				ARGU::ValueArg<float>	thArg(		"-t"			,"--threshold"				,"Threshold to print correlation [default=0.8]"				,"float"	, 0.8);
				argline.addArg(&thArg);
				/*ARGU::ValueArg<int>  mwArg( "-mw"				,"--max_window"	,"maximum extension of core length"		,"int"	,30);
				argline.addArg(&sfnArg);*/
				
				// Parse the command line
        argline.parseLine();
        
        // Get the environment variables
        score_matrix_path	= smpArg.getValue();
				if (score_matrix_path.empty()){
					cout << "Error : Environment variable missing : SCORE_MAT_PATH\n\n";
					argline.print_usage();
					exit(0);
				}
				
				// Get arguments
				mult_ali_fname			= mafArg.getValue();
				verbose							= verbArg.getValue();
				output_name         = outArg.getValue();
				statistic						= msArg.getValue();
        threshold		        = thArg.getValue();
				
        // Check the options which they are needed 
        if (mult_ali_fname.empty()){
					cout << "Error : A file name for multiple alignment is needed\n\n";
					argline.print_usage();
					exit(0);
				}
				if (output_name.empty()){
					int beg = mult_ali_fname.find_last_of('/');
					int end = mult_ali_fname.find_last_of('.');
				  output_name = mult_ali_fname.substr(beg+1, end - beg - 1);
				}
      } catch (exception &e) {
        throw;
      }
    }

public:
  
    /* List of options */
		string score_matrix_path;  	// The path to the scoring matrix
		string mult_ali_fname;			// The file name of the multiple alignment to read
		string output_name;         // Name of the output files
  	string statistic;           // Name of the statistic calculated by MstatX++
		bool verbose;								// Switch for verbose mode
  	float threshold;            // Threshold for correlation print
    /* Destructor */
    ~Options(){};
  
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
};

#endif