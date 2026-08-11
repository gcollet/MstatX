/* Copyright (c) 2012 Guillaume Collet
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
#include <stdexcept>

using namespace std;

#define MIN(x,y)  (x < y ? x : y)


/*
 * Calculate statistics based on the measure 
 * proposed by Capra and Singh (2007)
 * S = (1 - gapFreq) * (λ * R(p,r) + (1 - λ) * R(q,r))
 * with R(p,r) : \sum_{a=1}^{K} p(a) log (p(a)/r(a))
 * i.e. a Kullback-Leibler divergence with probability
 * p and background probability r
 */

/* Small constant added to a column's zero-probability symbols so no
 * log(0) is ever taken; subtracted back from the non-zero symbols so
 * the column's probabilities still sum to ~1 (see the pseudo_counts
 * computation below). Named here instead of repeating the same
 * pow(10.0,-6.0) three times. */
static const float PSEUDO_COUNT = 1e-6f;

void
JensenStat :: calculate(Msa & msa)
{
	/* Init size */
	string alphabet = msa.getAlphabet();
	int L = msa.getNcol();
	int N = msa.getNseq();
	int K = static_cast<int>(alphabet.size());
	
	/* Allocate proba array */
	std::vector<std::vector<float> > proba(L, std::vector<float>(K, 0.0f));

	/* Calculate Sequence Weights */
	const vector<float> & w = msa.getSeqWeights();

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
				proba[x][a] = PSEUDO_COUNT;
				nb_abs++;
			}
		}
		/* reduce by the pseudo counts in order to have sum-of-proba = 1 */
		float pseudo_counts = static_cast<float>(nb_abs) * PSEUDO_COUNT / static_cast<float>(K - nb_abs);
		for (int a(0); a < K; a++){
			if (proba[x][a] > PSEUDO_COUNT){
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
		col_stat.push_back((1 - (lambda * score_left + (1.0 - lambda) * score_right)) * (1 - (static_cast<float>(msa.getGap(x)) / static_cast<float>(N))));
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
}
