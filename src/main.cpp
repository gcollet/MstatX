#include <iostream>
#include <ctime>

#include "msa.h"
#include "options.h"
#include "statistic.h"
#include "scoring_matrix.h"

using namespace std;

int main (int argc, char **argv) {
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
	
	/* Read the multiple alignment */
	cout << "Read multiple alignment\n";
	Msa msa(Options::Get().mult_ali_fname);
	
	/* Calculate the statistics */
	Statistic * stat = StatisticFactory::CreateByName(Options::Get().statistic);
	stat->calculateStatistic(msa);
	
	delete stat;
	
	t2 = clock();		
	cout << "Mstatx computed in "<< (t2 - t1) / (double)CLOCKS_PER_SEC <<" seconds\n";
	return 0;
}
