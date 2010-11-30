/*
 *  wentropy.cpp
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 10/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "wentropy.h"
#include "options.h"

#include <cmath>
#include <fstream>
#include <algorithm>

#define MIN(x,y)  (x < y ? x : y)

WEntStat :: ~WEntStat()
{
	if (proba != NULL){
		for (int i(0); i < ncol; ++i) {
			if (proba[i] != NULL){
				free (proba[i]);	
			}
		}
		free(proba);
	}
}

/* Calculate the weight of sequence i in the msa 
 * by the formula from Henikoff & Henikoff (1994)
 * w_i = \frac{1}{L}\sum_{x=1}^{L}\frac{1}{k_x n_{x_i}}
 * We use these notations in the code below
 */
float 
WEntStat :: calcSeqWeight(Msa & msa, int i)
{
	int x, seq;
	int k;								    /**< number of symbol types in a column */
	int n;								    /**< number of occurence of aa[i][x] in column x */
	int L = msa.getNcol();    /**< number of columns (i.e. length of the alignment) */
	int nseq = msa.getNseq(); /**< number of rows (i.e. number of sequences in the alignment) */
	float w;						    	/**< weight of sequence i */
	
	w = 0.0;
	for	(x = 0; x < L; ++x){
		k = msa.getNtype(x);
		/* Calculate the number of aa[i][x] in current column */
		n = 0;
		for(seq = 0; seq < nseq; ++seq){
		  if (msa.getSymbol(i, x) == msa.getSymbol(seq, x)){
				n++;
			}
		}
		w += (float) 1 / (float) (n * k); 
	}
	w /= (float) L;
	return w ;
}

void
WEntStat :: calculateStatistic(Msa & msa)
{
	/* Init size */
	ncol = msa.getNcol();
	nseq = msa.getNseq();
	string alphabet = msa.getAlphabet();
	
	/* Allocate proba array */
	proba = (float **) calloc(ncol, sizeof(float*));
	if (proba == NULL){
		fprintf(stderr,"Cannot Allocate probability matrix\n");
		exit(0);
	}
	for (int i(0); i < ncol; i++){
		proba[i] = (float *) calloc(alphabet.size(), sizeof(float));
		if (proba[i] == NULL){
			fprintf(stderr,"Cannot Allocate probability submatrix\n");
			exit(0);
		}
	}
	
	/* Calculate Sequence Weights */
	for (int seq(0); seq < nseq; ++seq){
		seq_weight.push_back(calcSeqWeight(msa,seq));
	}
	
	/* Print if verbose mode on */
	if (Options::Get().verbose){
		cout << "Seq weights :\n";
		for (int seq(0); seq < nseq; ++seq){
		  cout.width(10);
		  cout << seq_weight[seq] << "\n";
		}
		cout << "\n";
	}
	
	/* Calculate aa proba and conservation score by columns */
	float lambda = 1.0 / log(MIN(alphabet.size(),nseq));
	
  for (int x(0); x < ncol; x++){
		col_cons.push_back(0.0);
		for (int a(0); a < alphabet.size(); a++){
		  for (int j(0); j < nseq; j++){
				if(msa.getSymbol(j, x) == alphabet[a]){
					proba[x][a] += seq_weight[j];
				}
			}
			if (proba[x][a] != 0.0){
				col_cons[x] -= proba[x][a] * log(proba[x][a]);
			}
		}
		col_cons[x] *= lambda;
	}
	
	cout << "Score is based on wentropy + gap counts\n";
	cout << "S = (1 - entropy) * (1 - gap_freq)\n";
	
	/* Print Conservation score in output file */
	ofstream file(Options::Get().output_name.c_str());
	if (!file.is_open()){
	  cerr << "Cannot open file " << Options::Get().output_name << "\n";
		exit(0);
	}
	for (int col(0); col < ncol; ++col){
	  file << (1.0 - col_cons[col]) * (1 - ((float) msa.getGap(col) / (float) nseq)) << "\n";
	}
	file.close();
}

