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

#include "mvector.h"
#include "options.h"
#include "scoring_matrix.h"

#include <cmath>
#include <fstream>

using namespace std;

float 
MVectStat :: normVect(vector<float> vect){
	float score= 0.0;
	for(int i(0); i < vect.size(); ++i){
		score += vect[i] * vect[i];
	}
	return sqrt(score);
}

void
MVectStat :: calculateStatistic(Msa & msa)
{
	ncol = msa.getNcol();
	nseq = msa.getNseq();
	/* Get the scoring matrix */
	ScoringMatrix score_mat(Options::Get().score_matrix_path + "/mclachlan71.mat");
	/* Remove the unknown symbol from msa (consider them as gaps)*/
	string sm_alphabet = score_mat.getAlphabet();
	msa.fitToAlphabet(sm_alphabet);
	
	/* Calculate the mean vector for each column */
	int alph_size = score_mat.getAlphabetSize();
	vector<vector<float> > means(ncol);
	vector<float> mean_col(alph_size, 0.0);
	for (int col(0); col < ncol; col++) {
		for (int seq(0); seq < nseq; ++seq) {
			if (msa.getSymbol(seq,col) == '-'){
				continue;
			} else {
				for (int a(0); a < alph_size; ++a) {
					mean_col[a] += score_mat.normScore(sm_alphabet[a],msa.getSymbol(seq,col));
				}
			}
		}
		for (int a(0); a < alph_size; ++a) {
			mean_col[a] /= (float) nseq;
		}
		means[col] = mean_col;
	}
	
	/* Print the output */
	ofstream file(Options::Get().output_name.c_str());
	if (!file.is_open()){
	  cerr << "Cannot open file " << Options::Get().output_name << "\n";
		exit(0);
	}
	for (int a(0); a < alph_size; ++a) {
		file << score_mat.getAlphabet()[a] << " ";
	}
	file << "\n";
	for (int col(0); col < ncol; col++) {
  	for (int a(0); a < alph_size; ++a) {
	  	mean_col[a] /= (float) nseq;
			file << means[col][a] << " ";
		}
		file << "\n";
	}
	
	
	
}