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
#include <iostream>
#include <fstream>
#include <cmath>

#include "options.h"
#include "scoring_matrix.h"

using namespace std;

/** Constructor from a filename fname.
 *  Matrices are all in format defined by AAindex web site :
 *  http://www.genome.jp/aaindex/
 */
ScoringMatrix :: ScoringMatrix(string fname)
{
  /* Open file */
	if(fname.empty()){
		cerr << "Error, score matrix file name is empty\n";
		exit(0);
	}
	if (Options::Get().verbose){
		cout << "Read Scoring Matrix in " << fname << "\n";
	}
	ifstream file(fname.c_str());
	if (!file.good()){
	  cerr << "Cannot open file " << fname << "\n";
		exit(0);
	}
	
	/* Read file */
	string s;
	getline(file,s);
	while (file.good() && s[0] != 'M'){
	  getline(file,s);
	}
	
	/* Read Alphabet */
	int alphabet_begin = s.find('=') + 2;
	int alphabet_size = s.find(',') - alphabet_begin;
	alphabet = s.substr(alphabet_begin,	alphabet_size);
	
	/* Allocate the matrix */
	matrix = (float **) calloc (alphabet_size, sizeof(float *));
	if (matrix == NULL){
	  cerr << "Cannot allocate scoring matrix\n";
		exit(0);
	}
	for (int i(0); i < alphabet_size; ++i) {
		matrix[i] = (float *) calloc(i + 1, sizeof(float));
		if (matrix[i] == NULL){
			cerr << "Cannot allocate scoring matrix[" << i << "]\n";
			exit(0);
		}
	}
	
	
	/* Read the matrix */
  min = 1000; max = -1000;
	for (int i(0); i < alphabet_size; ++i) {
		getline(file,s);
		for (int j(0); j <=i ; j++){
			matrix[i][j] = atof(s.substr(j*8, 8).c_str());
      if (matrix[i][j] < min)
        min = matrix[i][j];
      if (matrix[i][j] > max)
        max = matrix[i][j];
		}
	}
	
	
	/* Allocate the norm_matrix */
	norm_matrix = (float **) calloc (alphabet_size, sizeof(float *));
	if (norm_matrix == NULL){
	  cerr << "Cannot allocate the normalized amino acid vector\n";
		exit(0);
	}
	for (int i(0); i < alphabet_size; ++i) {
		norm_matrix[i] = (float *) calloc(i+1, sizeof(float));
		if (norm_matrix[i] == NULL){
			cerr << "Cannot allocate norm_matrix[" << i << "]\n";
			exit(0);
		}
	}
	
	/* Calculate normalized vector */
	
	for (int i(0); i < alphabet_size; ++i) {
		for (int j(0); j <= i; ++j){
			norm_matrix[i][j] = (score(alphabet[i], alphabet[j]) - min) / (max - min);
		}
	}
	
	if (Options::Get().verbose){
		cout << "Normalized :\n";
		for (int i(0); i < alphabet_size; ++i) {
			cout.width(9);
			cout << alphabet[i];
			for (int j(0); j <= i; j++){
				cout.width(9);
				cout << norm_matrix[i][j];
			}
			cout << "\n";
		}
		cout << "\n";
		cout.width(9);
		cout << ' ';
		for (int j(0); j <= alphabet_size ; j++){
			cout.width(9);
			cout << alphabet[j];
		} 
		cout << "\n\n";
	}
	
	is_set = true;
}

/* Destructor */
ScoringMatrix :: ~ScoringMatrix()
{
	if (is_set){
		for (int i(0); i < alphabet.size(); ++i){
			free(matrix[i]);
		}
		free(matrix);
		for (int i(0); i < alphabet.size(); ++i){
			free(norm_matrix[i]);
		}
		free(norm_matrix);
	}
	
}

int 
ScoringMatrix :: index(char aa)
{
	int pos = alphabet.find(aa);
	if (pos >= alphabet.size()){
		cerr << "Symbol " << aa << " is not in alphabet\n";
		exit(0);
	} 
	return pos;
}


float
ScoringMatrix :: score(char aa1, char aa2)
{
	int x,y;
  int pos1 = index(aa1);
	int pos2 = index(aa2);
	if (pos1 > pos2){
		x = pos1; y = pos2;
	} else {
	  x = pos2; y = pos1;	
	}
	return matrix[x][y];
}

float 
ScoringMatrix :: normScore(char aa1, char aa2)
{
	int x,y;
  int pos1 = index(aa1);
	int pos2 = index(aa2);
	if (pos1 > pos2){
		x = pos1; y = pos2;
	} else {
	  x = pos2; y = pos1;	
	}
	return norm_matrix[x][y];
}
