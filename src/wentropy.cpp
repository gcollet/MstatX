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

#include "wentropy.h"
#include "options.h"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm>

using namespace std;

#define MIN(x,y)  (x < y ? x : y)


/** calculate(Msa & msa)
 *
 * Calculate wentropy statistic and print it in the output file
 * The wentropy score is calculated as presented by Valdar (2002)
 * in equations (50), (51), and (52) :
 * For each column x : t(x) = \lambda_t \sum_{a \in K} p_a log(p_a)
 * With : \lambda_t = \frac{1.0}{log(min(K,N))}
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
	int K = static_cast<int>(alphabet.size());
	
	/* Allocate probabilities array */
	std::vector<std::vector<float> > p(L, std::vector<float>(K, 0.0f));
	
	/* Calculate Sequence Weights */
	const vector<float> & w = msa.getSeqWeights();
	
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
}
