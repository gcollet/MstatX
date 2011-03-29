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
				ARGU::ValueArg<string>	msArg(	"-s"			,"--statistic"				,"Name of the statistic used [default=trident]"						,"string"	,"wentropy");
				argline.addArg(&msArg);
				ARGU::ValueArg<float>		thArg(	"-t"			,"--threshold"				,"Threshold to print correlation [default=0.8]"				,"float"	, 0.8);
				argline.addArg(&thArg);
				ARGU::ValueArg<float>		taArg(	"-a"			,"--trident_a"				,"Factor applied to the first member of trident score  [default=1.0]"				,"float"	, 1.0);
				argline.addArg(&taArg);
				ARGU::ValueArg<float>		tbArg(	"-b"			,"--trident_b"				,"Factor applied to the second member of trident score [default=0.5]"				,"float"	, 0.5);
				argline.addArg(&tbArg);
				ARGU::ValueArg<float>		tcArg(	"-c"			,"--trident_c"				,"Factor applied to the third member of trident score  [default=3.0]"				,"float"	, 3.0);
				argline.addArg(&tcArg);
				ARGU::ValueArg<int>			nsArg(	"-n"			,"--nb_seq"						,"Maximum number of sequences read in the multiple alignment  [default=500]"				,"int"	, 500);
				argline.addArg(&nsArg);
				ARGU::ValueArg<int>			wArg(		"-w"			,"--window"						,"Number of side columns (jensen score)"		,"int"	,3);
				argline.addArg(&wArg);
				ARGU::SwitchArg					bArg(		"-ba"			,"--basic"					,"Add a basic output nbAA/col"																				,"string");
				argline.addArg(&bArg);
				
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
				factor_a						= taArg.getValue();
				factor_b						= tbArg.getValue();
				factor_c						= tcArg.getValue();
				nb_seq							= nsArg.getValue();
				window							=	wArg.getValue();
				basic								= bArg.getValue();
        // Check the options which they are needed 
        if (mult_ali_fname.empty()){
					cout << "Error : A file name for multiple alignment is needed\n\n";
					argline.print_usage();
					exit(0);
				}
				if (output_name.empty()){
					int beg = mult_ali_fname.find_last_of('/');
					int end = mult_ali_fname.find_last_of('.');
				  output_name = mult_ali_fname.substr(beg+1, end - beg - 1) + ".stat";
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
		bool basic;									// Switch for basic information output
  	float threshold;            // Threshold for correlation print
		float factor_a;							// Factor applied to the first  member of trident score
		float factor_b;             // Factor applied to the second member of trident score
		float factor_c;							// Factor applied to the third  member of trident score
		int nb_seq;
	  int window;                 // Size of the window to take in account side columns
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

