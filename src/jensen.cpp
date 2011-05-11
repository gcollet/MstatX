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

#include "jensen.h"
#include "options.h"
#include "scoring_matrix.h"

#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>

#define MIN(x,y)  (x < y ? x : y)


/* Calculate the weight of sequence i in the msa 
 * by the formula from Henikoff & Henikoff (1994)
 * w_i = \frac{1}{L}\sum_{x=1}^{L}\frac{1}{k_x n_{x_i}}
 * We use these notations in the code below
 */
float 
JensenStat :: calcSeqWeight(Msa & msa, int i)
{
	int x, seq;
	int k;								    /**< number of symbol types in a column */
	int n;								    /**< number of occurence of aa[i][x] in column x */
	int L = msa.getNcol();    /**< number of columns (i.e. length of the alignment) */
	int N = msa.getNseq(); /**< number of rows (i.e. number of sequences in the alignment) */
	float w;						    	/**< weight of sequence i */
	
	w = 0.0;
	for	(x = 0; x < L; ++x){
		k = msa.getNtype(x);
		/* Calculate the number of aa[i][x] in current column */
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

/*
 * Calculate statistics based on the measure 
 * proposed by Capra and Singh (2007)
 * S = (1 - gapFreq) * (λ * R(p,r) + (1 - λ) * R(q,r))
 * with R(p,r) : \sum_{a=1}^{K} p(a) log (p(a)/r(a))
 * i.e. a Kullback-Leibler divergence with probability
 * p and background probability r
 */
void
JensenStat :: calculate(Msa & msa)
{
	/* Init size */
	string alphabet = msa.getAlphabet();
	int L = msa.getNcol();
	int N = msa.getNseq();
	int K = alphabet.size();
	
	/* Allocate proba array */
	float **proba = (float **) calloc(L, sizeof(float*));
	if (proba == NULL){
		fprintf(stderr,"Cannot Allocate probability matrix\n");
		exit(0);
	}
	for (int i(0); i < L; i++){
		proba[i] = (float *) calloc(K, sizeof(float));
		if (proba[i] == NULL){
			fprintf(stderr,"Cannot Allocate probability submatrix\n");
			exit(0);
		}
	}

	/* Calculate Sequence Weights */
	vector<float> w;
	for (int seq(0); seq < N; ++seq){
		w.push_back(calcSeqWeight(msa,seq));
	}

	/* Background distribution of amino acids 
	 * These background frequencies are used in sca paper
	 */
	map<char,float> q;
	q['A'] = 0.073;
	q['C'] = 0.025;
	q['D'] = 0.050;
	q['E'] = 0.061;
	q['F'] = 0.042;
	q['G'] = 0.072;
	q['H'] = 0.023;
	q['I'] = 0.053;
	q['K'] = 0.064;
	q['L'] = 0.089;
	q['M'] = 0.023;
	q['N'] = 0.043;
	q['P'] = 0.052;
	q['Q'] = 0.040;
	q['R'] = 0.052;
	q['S'] = 0.073;
	q['T'] = 0.056;
	q['V'] = 0.063;
	q['W'] = 0.013;
	q['Y'] = 0.033;
	
	/* Calculate aa proba by columns */
	float lambda = 0.5;
	
  for (int x(0); x < L; ++x){
		int nb_abs = 0;
		for (int a(0); a < K; a++){
		  for (int j(0); j < N; ++j){
				if(msa.getSymbol(j, x) == alphabet[a]){
					proba[x][a] += w[j];
				}
			}
			if (proba[x][a] == 0.0){
				proba[x][a] = pow(10.0,-6.0);
				nb_abs++;
			}
		}
		/* reduce by the pseudo counts in order to have sum-of-proba = 1 */
		float pseudo_counts = (float) nb_abs * pow(10.0,-6.0) / (float) (K - nb_abs);
		for (int a(0); a < K; a++){
			if (proba[x][a] > pow(10.0,-6.0)){
				proba[x][a] -= pseudo_counts;
			}
		}
	}
	
	/* Calculate conservation score by columns */
	for (int x(0); x < L; ++x){
		float score_left = 0.0;
		float score_right = 0.0;
		for (int a(0); a < K; a++){
			char aa = alphabet[a];
			if (aa != '-' && aa != 'X' && aa != 'Z' && aa != 'B'){
				score_left  += proba[x][a] * log(proba[x][a] / (lambda * proba[x][a] + (1.0 - lambda) * q[aa]));
				score_right += q[aa] * log(q[aa] / (lambda * proba[x][a] + (1.0 - lambda) * q[aa]));
			}
		}
		col_stat.push_back((1 - (lambda * score_left + (1.0 - lambda) * score_right)) * (1 - ((float) msa.getGap(x) / (float) N)));
	}
	
	/* Add Side columns effect */
	/*int window = Options::Get().window;
	for (int x(0); x < L; ++x){
		float score = col_stat[x];
		float side_score = 0.0;
		for (int i(x-1); i >= 0 && i >= x - window; i--) {
			side_score += col_stat[i];
		}
		for (int i(x+1); i < L && i <= x + window; i++) {
			side_score += col_stat[i];
		}
		side_score /= (2 * window);
		cout << 0.5 * (score + side_score) <<"\n";
		col_stat[x] = 0.5 * (score + side_score);
	}*/
	
	cout << "\nScore is based on Jensen-Shannon measure\n";
	cout << "S = λ R(p,r) + (1 - λ) R(q,r)\n\n";
	
	
	for (int i(0); i < L; ++i) {
		free (proba[i]);	
	}
	free(proba);
}
