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
KabatStat :: calculate(Msa & msa)
{
	int k;	               /**< number of amino acid types in column */
	int n1;                /**< number of occurences of the most represented amino acid in a column */
	int N = msa.getNseq(); /**< number of sequences in the msa*/
	int ncol = msa.getNcol();  /**< number of columns in the multiple alignment */
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
		col_stat.push_back((float)k / (float) n1);
	}
	
	cout << "\nScore is based on kabat score\n";
	cout << "S = k / n1\n\n";
}
