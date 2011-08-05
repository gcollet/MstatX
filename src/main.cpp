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

#include <iostream>
#include <ctime>

#include "msa.h"
#include "options.h"
#include "statistic.h"
#include "scoring_matrix.h"

using namespace std;

int main (int argc, char **argv)
{
	clock_t t1,t2;
	t1 = clock();
	
	/* 
	 * Parses command line 
	 */
	try {
		Options::Get().Parse(argc, argv);
	}
	catch (exception &e) {
		cerr << e.what() << "\n";
		Options::Get().print_usage();
		exit(0);
	}
	
	/* 
	 * Initiates Statistic factory
	 */
	try {
	  AddAllStatistics();
	} catch (exception &e){ 
		cerr << e.what() << "\n";
		exit(0);
	}
	
	/*
	 * Read the multiple alignment 
	 */
	Msa msa(Options::Get().mult_ali_fname);
	
	/*
	 * If Basic, output basic information
	 */
	if (Options::Get().basic) {
		msa.printBasic();
	}
	
	/* 
	 * Calculate the statistic
	 */
	Statistic * stat = StatisticFactory::CreateByName(Options::Get().statistic);
	stat->calculate(msa);
	if (Options::Get().basic){
		msa.printBasic();
	} else {
		stat->print(msa);
	}
	delete stat;
	
	/*
	 * Print time
	 */
	t2 = clock();		
	cout << "Mstatx computed in "<< (t2 - t1) / (double)CLOCKS_PER_SEC <<" seconds\nResults are written in " << Options::Get().output_name << "\n\n";
	return 0;
}
