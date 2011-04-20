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
#include "mut_info.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#define MIN(x,y)  (x < y ? x : y)

using namespace std;

/** calcSeqWeight(Msa & msa, int i)
 *
 * Calculate the weight of sequence i in the msa 
 * by the formula from Henikoff & Henikoff (1994)
 * w_i = \frac{1}{L}\sum_{x=1}^{L}\frac{1}{k_x n_{x_i}} (LateX code)
 * We use these notations in the code below
 */
float 
MIStat :: calcSeqWeight(Msa & msa, int i)
{
	int x;										/**< used to parse msa columns */
	int seq;									/**< used to parse msa rows */
	int k;								    /**< number of symbol types in a column */
	int n;								    /**< number of occurence of aa[i][x] in column x */
	int L = msa.getNcol();    /**< number of columns (i.e. length of the alignment) */
	int N = msa.getNseq();    /**< number of rows (i.e. number of sequences in the alignment) */
	float w;						    	/**< weight of sequence i */
	
	w = 0.0;
	for	(x = 0; x < L; ++x){
		k = msa.getNtype(x);
		n = 0;
		for(seq = 0; seq < N; ++seq){
		  if (msa.getSymbol(i, x) == msa.getSymbol(seq, x)){
				n++;
			}
		}
		w += (float) 1 / (float) (n * k); 
	}
	w /= (float) L;
	return w ;
}


/** calculateStatistic(Msa & msa)
 *
 * Calculate Mutatual information statistic
 * First, a shannon-entropy score is calculated for each columns and for each pair of columns :
 * For each column x : t(x) = \lambda_t \sum_{a \in K} p_a log p_a
 * With : p_a = \sum_{i = 1}^{N}\left\{\begin{array}{l}w_i \mbox{ if }a=msa[i][x]\\0 \mbox{ else}\end{array}\right.
 *
 * Then the Mutual information of columns c and d is: MI_{cd} = t(c) + t(d) - t(c,d)
 *
 * These notations are used in the code
 */
void
MIStat :: calculateStatistic(Msa & msa)
{
	string alphabet = msa.getAlphabet();
	
	/* Init sizes */
	int L = msa.getNcol();
	int N = msa.getNseq();
	int K = alphabet.size();
	int npair = L * (L - 1) / 2;
	
	/* Allocate probabilities array */
	float ** p = (float **) calloc(L, sizeof(float*));
	if (p == NULL){
		fprintf(stderr,"Cannot Allocate probability matrix\n");
		exit(0);
	}
	for (int i(0); i < L; i++){
		p[i] = (float *) calloc(K, sizeof(float));
		if (p[i] == NULL){
			fprintf(stderr,"Cannot Allocate probability submatrix\n");
			exit(0);
		}
	}
	
	/* Calculate Sequence Weights */
	vector<float> w;
	for (int seq(0); seq < N; ++seq){
		w.push_back(calcSeqWeight(msa,seq));
	}
	
	/* Calculate aa proba and conservation score by columns */
	float lambda = 1.0 / log(MIN(K,N));
	float max = 0.0;
	vector<float> t;
	int x;
	for (x = 0; x < L; ++x){
		t.push_back(0.0);
		for (int a(0); a < K; ++a){
		  for (int j(0); j < N; ++j){
				if(msa.getSymbol(j, x) == alphabet[a]){
					p[x][a] += w[j];
				}
			}
			if (p[x][a] != 0.0){
				t[x] -= p[x][a] * log(p[x][a]);
			}
		}
		if (t[x] > max){
			max = t[x];
		}
		//t[x] *= lambda;
	}
	
	/* Calculate joint shannon-entropy values */
	vector< vector<float> > joint(L);
#pragma omp parallel for
	for  (x = 0; x < L-1; ++x) {
		joint[x] = vector<float>(L,0.0);
		for (int y(x+1); y < L; ++y) {
			for (int a(0); a < K; ++a) {
				for (int b(0); b < K; ++b) {
					float p = 0.0;
					for (int j(0); j < N; ++j) {
						if(msa.getSymbol(j, x) == alphabet[a] && msa.getSymbol(j, y) == alphabet[b]){
							p += w[j];
						}
					}
					if (p != 0.0){
						joint[x][y] -= p * log(p);
					}
				}
			}
			//joint[x][y] *= lambda * lambda;
		}
	}
	
	/* Calculate mutual information */
	mi = vector< vector<float> >(L);
	for  (int x(0); x < L-1; ++x) {
		mi[x] = vector<float>(L,0.0);
		for (int y(x+1); y < L; ++y) {
			mi[x][y] = (t[x] + t[y] - joint[x][y]) / max;
		}
	}
	
	/* Free all */
	for (int i(0); i < L; ++i) {
		free (p[i]);	
	}
	free(p);
}


/** printStatistic(Msa & msa)
 *
 * Print Mutual information scores in output file 
 */
void 
MIStat :: printStatistic(Msa & msa){
	ofstream file(Options::Get().output_name.c_str());
	if (!file.is_open()){
	  cerr << "Cannot open file " << Options::Get().output_name << "\n";
		exit(0);
	}
	for  (int x(0); x < mi.size() - 1; ++x) {
		for (int y(0); y < mi.size(); ++y) {
			if (y>x){
				file << mi[x][y] << " ";
			} else {
				file << 0.0 << " ";
			}
		}
		file << "\n";
	}
	file.close();
}


