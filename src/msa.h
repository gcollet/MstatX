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

#ifndef __MSA_H__
#define __MSA_H__

#include <vector>
#include <string>

using namespace std;

class Msa
{
protected:
	string alphabet;
	vector<string> mali_name;			/**< Name of sequences of the multiple alignment */
	vector<string> mali_seq;			/**< Sequences of the multiple alignment */
	vector<string> aa_type_list;	/**< List of aa type in each column (size = ncol * 20) */
	vector<int>    gap_counts;		/**< Number of gaps in each column */
	vector<float>  aa_freq;				/**< Frequency of amino acids types in the overall multiple alignment */
	vector<float>  entropy;				/**< Entropy of each column of the multiple alignment */
	vector<int>    nb_type;				/**< Number of amino acid types in the column */
	
	int nseq;											/**< Number of sequences in the multiple alignment */
	int ncol;											/**< Number of columns in the multiple alignment */
	
	void countGap();							/**< Count the number of gap in each column */
	void countFreq();							/**< Calculate the frequencies of each amino acid type in the multiple alignment */
	void countType();							/**< Calculate the number of different amino acid types in each column */
	void countEntropy();					/**< Calculate the entropy of each column in the multiple alignment */
	void defineAlphabet();				/**< Define the alphabet used in the multiple alignment */
	
public:
	Msa(string fname);
	~Msa(){};
	
	int   getAaPos(char aa);		/**< Converts a char in his position in alphabet */
	float getFreq(char aa);			/**< Return the frequency of amino acid aa in the overall multiple alignment */
	int   getGap(int col);			/**< Return the number of gaps in the column col */
	vector<int> getGapCount(){return gap_counts;};
	
	int   getNcol() const {return ncol;};									/**< Returns ncol value */
	int   getNseq() const {return nseq;};									/**< Returns nseq value */
	int   nbGap(int col) const {return gap_counts[col];};	/**< Returns the number of gaps in column col */
	bool  isInclude(string alph1);												/**< True if the alphabet of the multiple alignment is included in the alphabet alph1 */
	
	string getCol(int col);																/**< Returns a column as a string */
	string getAlphabet() const{return alphabet;};					/**< Returns the alphabet of the msa */
	
	char getSymbol(int seq, int col){return mali_seq[seq][col];};	/**< Return symbol row seq, column col */
	int getNtype(int col){return nb_type[col];};									/**< Return the number of different amino acids in the column col */
	string getTypeList(int col){return aa_type_list[col];};				/**< Return the list of amino acid types in the column col */
	
	void fitToAlphabet(string alph1);																		/**< if a symbol of the msa is not in alphabet alph1, then it is changed in a gap '-' */
	void printBasic();
};

#endif

