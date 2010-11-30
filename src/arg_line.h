/*
 *  arg_line.h
 *  ARGU
 *
 *  Created by Guillaume Collet on 21/08/10.
 *
 */

#ifndef _ARGLINE_H_
#define _ARGLINE_H_

#include "argument.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

using namespace std;

namespace ARGU {
	
	/** \class ArgLine
	 * \brief This class implements a command line parser
	 *
	 */
	class ArgLine 
	{
	private:
		string app_name;							/**< Name of the application */
		int argc;											/**< Size of argv */
		char * const * argv;					/**< array of strings in the command line */
		map<string, Arg *> arg_list;	/**< Map of argument objects by their flag */
		map<string, string> ref_list;	/**< Map of small flag / long flag correspondence */
	
		/** Reduce a pathname in a basename
		 */
		string basename(string fname){
			int pos = fname.find_last_of('/');
			return fname.substr(pos+1, fname.size() - pos);
		};
		
	public:
		/**
		 * \param _argc size of argv
		 * \param _argv array of strings in the command line
		 *
		 * The constructor creates a default argument -h/--help
		 */
		ArgLine(int _argc, char * const _argv[]):argc(_argc),argv(_argv)
		{
			addArg(new SwitchArg("-h", "--help", "Print this help", "string"));
		};
		
		/**
		 * The destructor deletes the default argument
		 */
		~ArgLine()
		{
		  delete arg_list["-h"];
		};
		
		/** 
		 * This function adds an argument to the arg_list
		 * and creates a correspondence between small
		 * and long flags
		 */
		void addArg(Arg * arg){
			arg_list[arg->getSmallFlag()] = arg;
			ref_list[arg->getSmallFlag()] = arg->getSmallFlag();
			ref_list[arg->getLongFlag()] = arg->getSmallFlag();
		}
		
		/**
		 * This function read the arguments in argv
		 * and sets the corresponding Arg objects
		 * in arg_list
		 */
		void parseLine(){
			string name;
			map<string,string>::iterator it;
			app_name = basename(argv[0]);
			for (int i(1); i < argc; ++i) { // for all strings in argv
				
				// find the correspondence in ref_list
				it = ref_list.find(argv[i]);
				if (it != ref_list.end()) {
					name = it->second;
				} else {
					cout << "argument " << argv[i] << " is not valid\n";
					exit(0);
				}
				
				// if help is activated then print usage
				if (name == "-h") {
					print_usage();
					exit(0);
				}
				
				// if a value is needed for this argument then read it
				if (arg_list[name]->needValue()){
					i++;
					arg_list[name]->setValue(argv[i]);
				} else {
					arg_list[name]->setValue(name);
				}
			}
		}
		
		/**
		 * This method prints the usage of the application
		 */
		void print_usage(){
			cout << "Usage : " << app_name << " -ma multiple_alignment\n\nOptions:\n";
			map<string,Arg *>::iterator it = arg_list.begin();
			map<string,Arg *>::iterator ite = arg_list.end();
			int sflag_size = 0;
			int lflag_size = 0;
			int desc_size = 0;
			while (it != arg_list.end()){
				if (it->second->getSmallFlag().size() > sflag_size) {
					sflag_size = it->second->getSmallFlag().size();
				}
				if (it->second->getLongFlag().size() > lflag_size) {
					lflag_size = it->second->getLongFlag().size();
				}
				if (it->second->getDescription().size() > desc_size) {
					desc_size = it->second->getDescription().size();
				}
				it++;
			}
			it = arg_list.begin();
			while (it != arg_list.end()){
				string flag = it->second->getSmallFlag() + ",";
				cout << "\t" << setw(sflag_size+1) << left << flag;
				cout << "\t\t" << setw(lflag_size) << left << it->second->getLongFlag()<< "\t\t: ";
				cout << setw(desc_size) << left << it->second->getDescription() << "\n";
				it++;
			}
		}
	};
}

#endif

