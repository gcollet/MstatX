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
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

#include "options.h"
#include "scoring_matrix.h"

/** Constructor from a filename fname.
 *  Matrices are all in format defined by AAindex web site :
 *  http://www.genome.jp/aaindex/
 */
ScoringMatrix :: ScoringMatrix(std::string fname)
{
  /* Open file */
	if(fname.empty()){
		throw std::runtime_error("score matrix file name is empty");
	}
	if (Options::Get().verbose){
		std::cout << "Read Scoring Matrix in " << fname << "\n";
	}
	std::ifstream file(fname.c_str());
	if (!file.good()){
		throw std::runtime_error("Cannot open file " + fname);
	}
	
	/* Read file */
	std::string s;
	getline(file,s);
	while (file.good() && s[0] != 'M'){
	  getline(file,s);
	}
	
	/* Read Alphabet */
	int alphabet_begin = (int) s.find('=') + 2;
	int alphabet_size = (int) s.find(',') - alphabet_begin;
	alphabet = s.substr(alphabet_begin,	alphabet_size);
	
	/* Allocate the matrix */
	matrix = (float **) calloc (alphabet_size, sizeof(float *));
	if (matrix == NULL){
		throw std::runtime_error("cannot allocate scoring matrix");
	}
	for (int i(0); i < alphabet_size; ++i) {
		matrix[i] = (float *) calloc(i + 1, sizeof(float));
		if (matrix[i] == NULL){
			throw std::runtime_error("cannot allocate scoring matrix");
		}
	}
	
	
	/* Read the matrix
	 * Each row is whitespace-separated (the AAindex format pads its
	 * columns to a fixed width, but that padding is not reliable enough
	 * to parse on: a value that prints one character shorter than its
	 * neighbours shifts every following fixed-width field on that row.
	 * Tokenizing on whitespace works regardless of column width, and
	 * also degrades gracefully on files from other matrix databases
	 * that don't pad at all. */
	min = 1000; max = -1000;
	for (int i(0); i < alphabet_size; ++i) {
		getline(file,s);
		std::istringstream row_stream(s);
		for (int j(0); j <=i ; j++){
			if (!(row_stream >> matrix[i][j])){
				throw std::runtime_error(
					"malformed scoring matrix row for symbol '" + std::string(1, alphabet[i]) +
					"': expected " + std::to_string(i + 1) + " values, could only read " + std::to_string(j));
			}
      if (matrix[i][j] < min)
        min = matrix[i][j];
      if (matrix[i][j] > max)
        max = matrix[i][j];
		}
	}
	
	
	/* Allocate the norm_matrix */
	norm_matrix = (float **) calloc (alphabet_size, sizeof(float *));
	if (norm_matrix == NULL){
		throw std::runtime_error("cannot allocate the normalized amino acid vector");
	}
	for (int i(0); i < alphabet_size; ++i) {
		norm_matrix[i] = (float *) calloc(i+1, sizeof(float));
		if (norm_matrix[i] == NULL){
			throw std::runtime_error("cannot allocate norm_matrix");
		}
	}
	
	/* Calculate normalized vector */
	
	for (int i(0); i < alphabet_size; ++i) {
		for (int j(0); j <= i; ++j){
			norm_matrix[i][j] = (score(alphabet[i], alphabet[j]) - min) / (max - min);
		}
	}
	
	if (Options::Get().verbose){
		std::cout << "Normalized :\n";
		for (int i(0); i < alphabet_size; ++i) {
			std::cout.width(9);
			std::cout << alphabet[i];
			for (int j(0); j <= i; j++){
				std::cout.width(9);
				std::cout << norm_matrix[i][j];
			}
			std::cout << "\n";
		}
		std::cout << "\n";
		std::cout.width(9);
		std::cout << ' ';
		for (int j(0); j <= alphabet_size ; j++){
			std::cout.width(9);
			std::cout << alphabet[j];
		} 
		std::cout << "\n\n";
	}
	
	is_set = true;
}

/* Destructor */
ScoringMatrix :: ~ScoringMatrix()
{
	if (is_set){
		for (int i(0); i < (int) alphabet.size(); ++i){
			free(matrix[i]);
		}
		free(matrix);
		for (int i(0); i < (int) alphabet.size(); ++i){
			free(norm_matrix[i]);
		}
		free(norm_matrix);
	}
	
}

int 
ScoringMatrix :: index(char aa)
{
	/* alphabet.find() returns std::string::npos when aa is absent, not
	 * some out-of-range int. Casting it to int (the previous check) wraps
	 * around to -1, which is never >= alphabet.size(): the absent-symbol
	 * error would silently never fire, and callers would go on to index
	 * the score matrix with -1 - undefined behavior. Comparing against
	 * std::string::npos directly, before any cast, avoids that. */
	std::string::size_type pos = alphabet.find(aa);
	if (pos == std::string::npos){
		throw std::runtime_error(std::string("symbol ") + aa + " is not in alphabet");
	} 
	return (int) pos;
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
