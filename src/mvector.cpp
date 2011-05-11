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
#include <iomanip>

using namespace std;

void
MVectStat :: calculate(Msa & msa)
{
	int L = msa.getNcol();
	int N = msa.getNseq();
	
	/* Get the scoring matrix */
	ScoringMatrix score_mat(Options::Get().score_matrix_path + "/" + Options::Get().score_matrix_fname);
	
	/* Remove the unknown symbol from msa (consider them as gaps)*/
	sm_alphabet = score_mat.getAlphabet();
	msa.fitToAlphabet(sm_alphabet);
	
	/* Calculate the mean vector for each column */
	int K = sm_alphabet.size();
	means = vector<vector<float> >(L);
	
	for (int col(0); col < L; col++) {
		vector<float> mean_col(K, 0.0);
		for (int seq(0); seq < N; ++seq) {
			if (msa.getSymbol(seq,col) == '-'){
				continue;
			} else {
				for (int a(0); a < K; ++a) {
					mean_col[a] += score_mat.normScore(sm_alphabet[a],msa.getSymbol(seq,col));
				}
			}
		}
		for (int a(0); a < K; ++a) {
			mean_col[a] /= (float) N;
		}
		means[col] = mean_col;
	}
}

void
MVectStat :: print(Msa & msa)
{
	/* Print the output */
	ofstream file(Options::Get().output_name.c_str());
	if (!file.is_open()){
	  cerr << "Cannot open file " << Options::Get().output_name << "\n";
		exit(0);
	}
	int K = sm_alphabet.size();
	file.precision(3);
	file << setw(10) << " ";
	for (int a(0); a < K; ++a) {
		file.precision(3);
		file << setw(10) << sm_alphabet[a];
	}
	file << "\n";
	for (int col(0); col < means.size(); col++) {
		file.precision(3);
		file << setw(10) << col + 1;
  	for (int a(0); a < K; ++a) {
			file.precision(3);	
			file << setw(10) << means[col][a];
		}
		file << "\n";
	}
	file.close();
}
