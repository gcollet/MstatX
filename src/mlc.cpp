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

#include <cmath>
#include <vector>
#include <fstream>

//#include "omp.h"
#include "mlc.h"
#include "msa.h"
#include "options.h"
#include "scoring_matrix.h"

MlcStat :: ~MlcStat(){
	if (pair_mat != NULL){
	  for (int i(0); i < ncol; ++i){
		  if (pair_mat[i] != NULL){
			  free(pair_mat[i]);	
			}
		}
		free(pair_mat);
	}
	if (cor_mat != NULL){
	  for (int i(0); i < mean.size(); ++i){
		  if (cor_mat[i] != NULL){
			  free(cor_mat[i]);	
			}
		}
		free(cor_mat);
	}
}

void 
MlcStat :: calculate(Msa & msa)
{
  /* Create the McLachlan score matrix */
	ScoringMatrix score_mat(Options::Get().score_matrix_path + "/" + Options::Get().score_matrix_fname);
	
	/* If alphabets are not matching, stop */
	if (!msa.isInclude(score_mat.getAlphabet())){
		cerr << "Error : Multiple alignment alphabet is not included in score matrix alphabet\n";
		cerr << "Score matrix alphabet       : " << score_mat.getAlphabet() << "\n";
		cerr << "Multiple alignment alphabet : " << msa.getAlphabet() << "\n";
		exit(0);
	}
	
	/* Initialization */
	int index = 0;
	ncol = msa.getNcol();
	int nseq = msa.getNseq();
	int max_gap = nseq / 10;
	vector<int> gap_count = msa.getGapCount();
	mean  = vector<float>(ncol,0.0);
	stdev = vector<float>(ncol,0.0);
	
	/* Alloc pair_mat */
	pair_mat = (float **) calloc (ncol, sizeof(float*));
	for (int i(0); i < ncol; ++i){
		if (msa.nbGap(i) < max_gap){
			pair_mat[i] = (float*) calloc (nseq * nseq, sizeof(float));
		} else {
		  pair_mat[i] = NULL;	
		}
	}
	
	/* Fill pair_mat + calculate mean and standard deviation */
	cerr<<"  Fill pair matrix \n";
	int col;

	for (col = 0; col < ncol; ++col){
		
		/* If the number of gap is not too high */
		if (msa.nbGap(col) < max_gap){
			
			/* Fill the column pair matrix */
			string column = msa.getCol(col);
			for (int k(0); k < nseq; ++k) {
	  		for (int l(0); l < nseq; ++l) {
				  index = k * nseq + l;
					if (column[k] != ' ' && column[k] != '-' && column[l] != ' ' && column[l] != '-' ) {
						pair_mat[col][index] = score_mat.score(column[k], column[l]);
					} else {
						pair_mat[col][index] = 0.0;
					}
				}
			}
			
			/* Calculate the mean of the column */
			for (int i(0); i < nseq * nseq; ++i){
				mean[col] += pair_mat[col][i];
			}
			mean[col] /= (float) ((nseq - msa.nbGap(col)) * (nseq - msa.nbGap(col)));
		
			/* Calculate the standard deviation of the column */
			float total = 0.0;
			for (int i(0); i < nseq * nseq; ++i){
			  total += (pair_mat[col][i] - mean[col]) * (pair_mat[col][i] - mean[col]);
			}
			stdev[col] = sqrt (total / (float) ((nseq - msa.nbGap(col)) * (nseq - msa.nbGap(col))));
		}
	}
	
	/* Allocate cor_mat */
	cor_mat = (float **) calloc(ncol, sizeof(float*));
	for (int col(0); col < ncol; ++col){
		cor_mat[col] = (float*) calloc(col + 1,sizeof(float));
	}
	
	vector<string> columns;
	for (int col1 = 0; col1 < ncol; ++col1){
	  columns.push_back(msa.getCol(col1));	
	}
	
	/* Calculate MLC coefficients */
	cout << "  Calculate MLC correlations\n";
	int col1;
//#pragma omp parallel for schedule (dynamic)
	for (col1 = 0; col1 < ncol; ++col1){
		if (stdev[col1] == 0.0 || gap_count[col1] > max_gap){
			continue;
		}
		int col2;
	  for (col2 = 0; col2 <= col1; ++col2){
			if (col1 == col2 || gap_count[col2] > max_gap || stdev[col2] == 0.0){
				continue;
			}
			/* Count non identical positions */
			int non_id = 0;             
      for (int seq(0); seq < nseq; seq++){ 
				if ( columns[col1][seq] != columns[col2][seq]) 
					non_id++;
      }
			/* If two columns are identical then, the score is 1.00 */
      if (non_id == 0) {
				cor_mat[col1][col2] = 1.00;
				continue;
      }
			int k;
			cor_mat[col1][col2] = 0.0;
			for (k = 0; k < nseq * nseq; k++) {
				cor_mat[col1][col2] += (pair_mat[col1][k] - mean[col1]) * (pair_mat[col2][k] - mean[col2]);
      }
      cor_mat[col1][col2] /= stdev[col1] * stdev[col2] * (float) (nseq * nseq) ;
		}
	}

}

void 
MlcStat :: print(Msa & msa)
{
	/* Print Conservation score in output file */
	ofstream file(Options::Get().output_name.c_str());
	if (!file.is_open()){
	  cerr << "Cannot open file " << Options::Get().output_name << "\n";
		exit(0);
	}
	
	for (int col1(0); col1 < ncol; ++col1){
		for (int col2(0); col2 <= col1; ++col2){
			if (cor_mat[col1][col2] > Options::Get().threshold){
				file << "pair ("<<col1<<","<<col2<<") = "<<cor_mat[col1][col2]<<" above "<<Options::Get().threshold<<"\n";	
			}
		}
	}
	file.close();
}
