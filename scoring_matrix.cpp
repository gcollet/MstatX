/*
 *  scoring_matrix.cpp
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 04/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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
	for (int i(0); i < alphabet_size; ++i) {
		getline(file,s);
		for (int j(0); j <=i ; j++){
			matrix[i][j] = atof(s.substr(j*8, 8).c_str());
		}
	}
	
	if (Options::Get().verbose){
		cout << "Alphabet : " << alphabet << " (size = " << alphabet_size << ")\n\n";
		cout << "Matrix :\n";
		for (int i(0); i < alphabet_size; ++i) {
			cout.width(9);
			cout << alphabet[i];
			for (int j(0); j <=i ; j++){
				cout.width(9);
				cout << matrix[i][j];
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
	
	/* Allocate the aa_vect */
	aa_vect = (float **) calloc (alphabet_size, sizeof(float *));
	if (aa_vect == NULL){
	  cerr << "Cannot allocate the normalized amino acid vector\n";
		exit(0);
	}
	for (int i(0); i < alphabet_size; ++i) {
		aa_vect[i] = (float *) calloc(i+1, sizeof(float));
		if (aa_vect[i] == NULL){
			cerr << "Cannot allocate aa_vect[" << i << "]\n";
			exit(0);
		}
	}
	
	/* Calculate normalized vector */
	for (int i(0); i < alphabet_size; ++i) {
		for (int j(0); j <= i; ++j){
			aa_vect[i][j] = score(alphabet[i], alphabet[j]) / sqrt(score(alphabet[i], alphabet[i]) * score(alphabet[j], alphabet[j]));
		}
	}
	
	if (Options::Get().verbose){
		cout << "Normalized :\n";
		for (int i(0); i < alphabet_size; ++i) {
			cout.width(9);
			cout << alphabet[i];
			for (int j(0); j <= i; j++){
				cout.width(9);
				cout << aa_vect[i][j];
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
			free(aa_vect[i]);
		}
		free(aa_vect);
	}
	
}

int 
ScoringMatrix :: index(char aa)
{
	int pos = alphabet.find(aa);
	if (pos >= alphabet.size()){
		cerr << "Symbol '"<< aa << "' not in alphabet\n";
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
