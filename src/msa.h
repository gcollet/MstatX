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

#pragma once

#include <array>
#include <vector>
#include <string>

class Msa
{
protected:
	std::string alphabet;
	std::array<int,256> alpha_index;	/**< alpha_index[static_cast<unsigned char>(c)] = position of c in `alphabet`, or -1. O(1) replacement for alphabet.find(c) */
	std::vector<std::string> mali_name;			/**< Name of sequences of the multiple alignment */
	std::vector<std::string> mali_seq;			/**< Sequences of the multiple alignment */
	std::vector<std::string> aa_type_list;	/**< List of aa type in each column (size = ncol * 20) */
	std::vector<int>    gap_counts;		/**< Number of gaps in each column */
	std::vector<float>  aa_freq;				/**< Frequency of amino acids types in the overall multiple alignment */
	std::vector<float>  entropy;				/**< Entropy of each column of the multiple alignment */
	std::vector<int>    nb_type;				/**< Number of amino acid types in the column */
	std::vector<float>  seq_weight;		/**< Cache for the Henikoff & Henikoff sequence weights, see getSeqWeights() */
	bool           seq_weight_computed;
	
	int nseq;											/**< Number of sequences in the multiple alignment */
	int ncol;											/**< Number of columns in the multiple alignment */
	
	void countGap();							/**< Count the number of gap in each column */
	void countFreq();							/**< Calculate the frequencies of each amino acid type in the multiple alignment */
	void countType();							/**< Calculate the number of different amino acid types in each column */
	void countEntropy();					/**< Calculate the entropy of each column in the multiple alignment */
	void defineAlphabet();				/**< Define the alphabet used in the multiple alignment */
	void rebuildAlphaIndex();		/**< Rebuild alpha_index to match the current `alphabet` string */
	
public:
	explicit Msa(const std::string & fname);
	~Msa() = default;
	
	int   getAaPos(char aa);		/**< Converts a char in his position in alphabet */
	float getFreq(char aa);			/**< Return the frequency of amino acid aa in the overall multiple alignment */
	int   getGap(int col);			/**< Return the number of gaps in the column col */
	std::vector<int> getGapCount() const {return gap_counts;};
	
	int   getNcol() const {return ncol;};									/**< Returns ncol value */
	int   getNseq() const {return nseq;};									/**< Returns nseq value */
	int   nbGap(int col) const {return gap_counts[col];};	/**< Returns the number of gaps in column col */
	bool  isInclude(const std::string & alph1);												/**< True if the alphabet of the multiple alignment is included in the alphabet alph1 */
	
	std::string getCol(int col);																/**< Returns a column as a string */
	std::string getAlphabet() const{return alphabet;};					/**< Returns the alphabet of the msa */
	
	char getSymbol(int seq, int col){return mali_seq[seq][col];};	/**< Return symbol row seq, column col */
	int getNtype(int col){return nb_type[col];};									/**< Return the number of different amino acids in the column col */
	std::string getTypeList(int col){return aa_type_list[col];};				/**< Return the list of amino acid types in the column col */
	
	void fitToAlphabet(const std::string & alph1);																		/**< if a symbol of the msa is not in alphabet alph1, then it is changed in a gap '-' */
	void printBasic();
	
	const std::vector<float> & getSeqWeights();		/**< Henikoff & Henikoff (1994) sequence weights, computed once in O(nseq*ncol) and cached */
};

