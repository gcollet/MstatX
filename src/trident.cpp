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

#include "trident.h"
#include "options.h"
#include "scoring_matrix.h"

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
TridStat :: calcSeqWeight(Msa & msa, int i)
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
  //cerr << "w = " << w << "\n";
	return w ;
}

/** normVect(vector<float> vect)
 *
 * Return the vector norm = √(∑v*v)
 */
float 
TridStat :: normVect(vector<float> vect){
	float score= 0.0;
	for(int i(0); i < vect.size(); ++i){
		score += vect[i] * vect[i];
	}
	return sqrt(score);
}



void
TridStat :: calculate(Msa & msa)
{
	/* Declare the vectors */
	vector<float> w;					/**< Weight of each sequence in the msa (size = nb sequences) */
	vector<float> t;					/**< t(x) = Shannon entropy score  + Weighted sequence Score */
	vector<float> r;					/**< r(x) = Stereochemical score */
	vector<float> g;					/**< g(x) = Gap Score */
	
	/* Init size */
	int L = msa.getNcol();
	int N = msa.getNseq();
	string alphabet = msa.getAlphabet();
	int K = alphabet.size();
	
  //cerr << "Seq Weights\n";
	/* Calculate Sequence Weights */
	for (int seq(0); seq < N; ++seq){
		w.push_back(calcSeqWeight(msa,seq));
	}

	/* Calculate t(x) = \frac{\sum_{a=1}^{K}p_a log(p_a)}{log(min(N,K))}
	 *						p_a = \sum_{i \in \{i|s(i) = a\}} w_i
	 *						w_i = \frac{1}{L} \sum_{x=1}^{L}\frac{1}{K_x n_{x_i}}
	 * Like in wentropy
	 */
	float lambda = 1.0 / log(MIN(K,N));
	
  for (int x(0); x < L; x++){
		t.push_back(0.0);
		for (int a(0); a < K; a++){
			float tmp_proba(0.0);
		  for (int j(0); j < N; j++){
				if(msa.getSymbol(j, x) == alphabet[a]){
					tmp_proba += w[j];
				}
			}
			if (tmp_proba != 0.0){
				t[x] -= tmp_proba * log(tmp_proba);
			}
		}
		t[x] *= lambda;
    //cerr << "t[" << x << "] = " << t[x] << "\n";
	}
	
  
	/* Calculate r(x) = \lambda_r \frac{1}{k_x}\sum_{a=1}^{k_x}|\bar{X}(x) - X_a|
	 *      \lambda_r = \frac{1}{\sqrt{20(max(M)-min(M))^2}}
	 *					  X_a = \left[ \begin{array}{c}M(a,a_1)\\M(a,a_2)\\.\\.\\.\\M(a,a_{20})\end{array}\right]
	 *							M is a normalized scoring matrix
	 */
	ScoringMatrix score_mat(Options::Get().score_matrix_path + "/" + Options::Get().score_matrix_fname);
	int alph_size = score_mat.getAlphabetSize();
	string sm_alphabet = score_mat.getAlphabet();
	
	//msa.fitToAlphabet(sm_alphabet);
	
	for (int x(0); x < L; x++){
		
		int ntype = msa.getNtype(x);
		string type_list = msa.getTypeList(x);
		
		int pos = type_list.find('-');
		if (pos < type_list.size()){
			type_list.erase(type_list.begin()+pos);
			ntype--;
		}
		if (ntype){
			/* Calculate Mean vector */
			vector<float> mean(alph_size, 0.0);
			for (int i(0); i < ntype; ++i){
				for (int a(0); a < alph_size; ++a){
					mean[a] += score_mat.normScore(sm_alphabet[a], type_list[i]);
				}
			}
			for (int a(0); a < alph_size; ++a){
				mean[a] /= ntype;
			}
		
			/* Calculate Score */
			float lambda = sqrt(alph_size * (score_mat.getMax() - score_mat.getMin()) * (score_mat.getMax() - score_mat.getMin()));
			float tmp_score = 0.0;
			for (int i(0); i < ntype; ++i){
				vector<float> diff_vect;
				for(int a(0); a < alph_size; ++a){
					diff_vect.push_back(mean[a] - score_mat.normScore(sm_alphabet[a], type_list[i]));
				}
				tmp_score += normVect(diff_vect);
			}
			tmp_score /= ntype;
			tmp_score /= lambda;
			r.push_back(tmp_score);
		} else {
			r.push_back(1.0);
		}
    //cerr << "r[" << x << "] = " << r[x] << "\n";
	}

	/* Calculate g(x) = nb_gap / nb_seq 
	 * Represents the proportion of gaps in the column
	 */
	for (int x(0); x < L; x++){
		g.push_back((float) msa.getGap(x) / (float) N);
    //cerr << "g[" << x << "] = " << g[x] << "\n";
	}
	
	for (int x(0); x < L; x++){
		col_stat.push_back(pow((1-t[x]),Options::Get().factor_a)*pow((1-r[x]),Options::Get().factor_b)*pow((1-g[x]),Options::Get().factor_c));
	}
	
}
