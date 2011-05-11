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

#include "wentropy.h"
#include "options.h"

#include <cmath>
#include <fstream>
#include <algorithm>

#define MIN(x,y)  (x < y ? x : y)


/** calcSeqWeight(Msa & msa, int i)
 *
 * Calculate the weight of sequence i in the msa 
 * by the formula from Henikoff & Henikoff (1994)
 * w_i = \frac{1}{L}\sum_{x=1}^{L}\frac{1}{k_x n_{x_i}} (LateX code)
 * We use these notations in the code below
 */
float 
WEntStat :: calcSeqWeight(Msa & msa, int i)
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


/** calculate(Msa & msa)
 *
 * Calculate wentropy statistic and print it in the output file
 * The wentropy score is calculated as presented by Valdar (2002)
 * in equations (50), (51), and (52) :
 * For each column x : t(x) = \lambda_t \sum_{a \in K} p_a log p_a
 * With : p_a = \sum_{i = 1}^{N}\left\{\begin{array}{l}w_i \mbox{ if }a=msa[i][x]\\0 \mbox{ else}\end{array}\right.
 *
 * These notations are used in the code
 */
void
WEntStat :: calculate(Msa & msa)
{
	string alphabet = msa.getAlphabet();
	
	/* Init sizes */
	int L = msa.getNcol();
	int N = msa.getNseq();
	int K = alphabet.size();
	
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
	
  for (int x(0); x < L; ++x){
		col_stat.push_back(0.0);
		for (int a(0); a < K; ++a){
		  for (int j(0); j < N; ++j){
				if(msa.getSymbol(j, x) == alphabet[a]){
					p[x][a] += w[j];
				}
			}
			if (p[x][a] != 0.0){
				col_stat[x] -= p[x][a] * log(p[x][a]);
			}
		}
		col_stat[x] *= lambda;
	}
	
	for (int i(0); i < L; ++i) {
		free (p[i]);	
	}
	free(p);
}
