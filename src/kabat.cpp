/*
 *  kabat.cpp
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 16/12/10.
 *
 */

#include "options.h"
#include "kabat.h"

#include <cmath>
#include <fstream>

/* Calculate the conservation of columns in the msa 
 * by the formula from Wu & Kabat (1970)
 * V(i) = k / n1 (in the original function, this score is
 * multiplicated by N)
 * We use these notations in the code below
 */
void
KabatStat :: calculateStatistic(Msa & msa)
{
	int k;	/**< number of amino acid types in column */
	int n1; /**< number of occurences of the most represented amino acid in a column */
	int N = msa.getNseq();	/**< number of sequences in the msa*/
	int ncol = msa.getNcol(); /**< number of columns in the multiple alignment */
	vector<int> naa;
	
	for	(int x(0); x < ncol; ++x){
		
		k = msa.getNtype(x);
		
		/* Find the most represented amino acid type (n1) */
		string aa_type = msa.getTypeList(x);
		vector<int> nb_aa(aa_type.size(), 0);
		for (int s(0); s < N; ++s){
			nb_aa[aa_type.find(msa.getSymbol(s, x))]++;
		}
		n1 = 0;
		for (int i(0); i < nb_aa.size(); ++i){
			if (nb_aa[i] > n1)
				n1 = nb_aa[i];
		}
		
		/* Calculate conservation from Wu & Kabat formula */
		col_cons.push_back((float)k / (float) n1);
	}
	
	cout << "\nScore is based on kabat score\n";
	cout << "S = k / n1\n\n";
	
	/* Print Conservation score in output file */
	ofstream file(Options::Get().output_name.c_str());
	if (!file.is_open()){
	  cerr << "Cannot open file " << Options::Get().output_name << "\n";
		exit(0);
	}
	for (int col(0); col < ncol; ++col){
	  file << col_cons[col] << "\n";
	}
	file.close();
}