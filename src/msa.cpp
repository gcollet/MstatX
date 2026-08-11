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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

#include "msa.h"
#include "options.h"

using namespace std;


/**************************************************************
 * This constructor of a multiple alignment reads the 
 * multiple alignment in a multi-fasta format.
 * Once read, the multiple alignment is analysed to find
 * the alphabet used, the number of gaps and the entropy of 
 * each column, and the frequency of each amino acid type.
 **************************************************************/
Msa :: Msa(const std::string & fname)
{
	seq_weight_computed = false;
	alpha_index.fill(-1);
	
	/* Open file */
	if (Options::Get().verbose){
		std::cout << "Read Multiple Alignment in " << fname << "\n";
	}
	std::ifstream file(fname.c_str());
	if (!file.good()){
		throw std::runtime_error("Cannot open file " + fname);
	}

	
	/* Read file */
	std::string s, tmp_seq;
	while (file.good() && static_cast<int>(mali_seq.size()) < Options::Get().nb_seq){
		getline(file,s);
		if (s[0] == '>'){
			if (mali_name.size() != 0){
			  mali_seq.push_back(tmp_seq);	
			}
			if (static_cast<int>(mali_seq.size()) < Options::Get().nb_seq){
  			mali_name.push_back(s.substr(1, s.find_first_of(' ') - 1));
			}
			tmp_seq.clear();
		} else {
			tmp_seq = tmp_seq + s;
		}
	}
	if (static_cast<int>(mali_seq.size()) < Options::Get().nb_seq){
	  mali_seq.push_back(tmp_seq);
	}
	
	nseq = static_cast<int>(mali_name.size());
	ncol = static_cast<int>(mali_seq[0].size());
	std::cout << "\nMultiple alignment : nb seq = "<<nseq<<", nb col = "<<ncol<<"\n";
	
	/* Change all mali.seq in upper case*/
	for (int i(0); i < nseq; ++i){
		for(int j(0); j < ncol; ++j){
			mali_seq[i][j] = toupper(mali_seq[i][j]);
		}
	}
	
	/* Analyse the multiple alignment */
	defineAlphabet();
	countGap();
	countFreq();
	countType();
	countEntropy();
	
	
	/* Print if verbose mode */
	if (Options::Get().verbose){
		cout << "\nAlphabet :\n";
		for (char c : alphabet){
			cout << c << ";";
		}
		cout << "\n";
		cout << "\nMultiple Alignment :\n";
		for (const auto & seq : mali_seq){
			cout << seq << "\n";
		}
		cout << "\nAA Frequencies :\n";
		for (float f : aa_freq){
			cout << f << ";";
		}
		cout << "\n";
		cout << "\nGap counts :\n";
		for (int g : gap_counts){
			cout << g << ";";
		}
		cout << "\n";
		cout << "\nAA Entropy :\n";
		for (int i(0); i < static_cast<int>(entropy.size()); ++i){
			if (gap_counts[i] < nseq/10){
		  	cout << entropy[i] << ";";
			} else {
				cout << "-12.0;";
			}
		}
		cout << "\n";
		cout << "\nAA Types :\n";
		for (int n : nb_type){
				cout << n << ";";
		}
		cout << "\n";
	}
}


/**************************************************************
 * countGap() calculate the number of gaps in each column
 **************************************************************/
void
Msa :: countGap(){
	for(int col(0); col < ncol; ++col){
		int gap = 0;
		for(int row(0); row < nseq; ++row){
			if (mali_seq[row][col] == '-' || mali_seq[row][col] == ' '){
				gap++;
			}
		}
		gap_counts.push_back(gap);
	}
}


/**************************************************************
 * countFreq() calculates the frequency of each amino acid
 * type in the overall multiple alignment
 **************************************************************/
void
Msa :: countFreq(){
	int total = 0;
	std::vector<int> tmp_freq(alphabet.size(), 0);
	
	aa_freq = std::vector<float>(alphabet.size());
	/* Count the number of each amino acid type defined in alphabet */
	for(int col(0); col < ncol; ++col){
		for(int row(0); row < nseq; ++row){
			if (mali_seq[row][col] != '-' && mali_seq[row][col] != ' '){
				total++;
			}
			int pos = alpha_index[static_cast<unsigned char>(mali_seq[row][col])];
			if (pos < 0){
				throw std::runtime_error("symbol is not in the alphabet");
			}
			tmp_freq[pos]++;
		}
	}

	/* Divide by the total */
	for (int i(0); i < static_cast<int>(aa_freq.size()); ++i){
		aa_freq[i] = static_cast<float>(tmp_freq[i]) ;
		aa_freq[i] /= static_cast<float>(total);
	}
}

/**************************************************************
 * countType() calculates the number of different types in
 * each column of the alignment
 **************************************************************/
void
Msa :: countType(){
	std::string aa_types;
	for(int col(0); col < ncol; ++col){
		aa_types.clear();
		for(int row(0); row < nseq; ++row){
			if (aa_types.find(mali_seq[row][col]) >= aa_types.size()){
			  aa_types.push_back(mali_seq[row][col]);
			}
		}
		aa_type_list.push_back(aa_types);
		nb_type.push_back(static_cast<int>(aa_types.size()));
	}
}

/**************************************************************
 * defineAlphabet() reads the multiple alignment to
 * determine all the symbols used in
 **************************************************************/
void
Msa :: defineAlphabet(){
	alphabet.clear();
	array<bool,256> seen;
	seen.fill(false);
	for(int col(0); col < ncol; ++col){
		for(int row(0); row < nseq; ++row){
			unsigned char c = mali_seq[row][col];
			if (!seen[c]){
				seen[c] = true;
				alphabet.push_back(static_cast<char>(c));
			}
		}
	}
	rebuildAlphaIndex();
}

/**************************************************************
 * rebuildAlphaIndex() rebuilds the O(1) char -> position
 * lookup table to match the current content of `alphabet`.
 * Must be called every time `alphabet` is mutated.
 **************************************************************/
void
Msa :: rebuildAlphaIndex(){
	alpha_index.fill(-1);
	for (int i(0); i < static_cast<int>(alphabet.size()); ++i){
		alpha_index[static_cast<unsigned char>(alphabet[i])] = i;
	}
}


/**************************************************************
 * getFreq(aa) returns the frequency of amino acid aa
 * in the overall multiple alignment
 **************************************************************/
float 
Msa :: getFreq(char aa) const {
  int pos = getAaPos(aa);
	if (pos == -1){
		throw std::runtime_error("symbol not in alphabet");
	}
	return aa_freq[pos];
}


/**************************************************************
 * getAaPos(aa) returns the index of aa in alphabet
 **************************************************************/
int
Msa :: getAaPos(char aa) const {
	return alpha_index[static_cast<unsigned char>(aa)];
};


/**************************************************************
 * getGap(col) returns the number of gap in column col
 **************************************************************/
int
Msa :: getGap(int col){
	return gap_counts[col];
}


/**************************************************************
 * countEntropy calculates the entropy of each column
 * by the following formula (Normalized Shannon Entropy):
 * entropy[col] = - \frac{1}{log(K)} \sum_{a=1}^{K}p_a log(p_a)
 * with 
 * K = alphabet length
 * p_a = probability to see amino acid of type a in the column
 * p_a = frequency of amino acid a in the column (nb_a / nseq)
 **************************************************************/
void 
Msa :: countEntropy(){
	entropy = std::vector<float>(ncol,0.0);
 
  for(int col(0); col < ncol; ++col){
		std::vector<float> lfreq(alphabet.size(), 0.0);
		for(int row(0); row < nseq; ++row){
			lfreq[getAaPos(mali_seq[row][col])] += 1.0;
		}
		for (float & f : lfreq){
		  f /= static_cast<float>(nseq);
			if (f > 0.0){
				if (f == 1.0){
				  entropy[col] = 0.0;	
				} else {
				  entropy[col] -= f * log(f);	
				}
			}
		}
		entropy[col] /= log(alphabet.size()-1); /* -1 because gaps are in the alphabet */
	}
}


/**************************************************************
 * isInclude(alph1) returns true if the alphabet of the
 * multiple alignment is include in the alphabet alph1
 **************************************************************/
bool
Msa :: isInclude(const std::string & alph1) const {
  for (char c : alphabet){
		if (alph1.find(c) == std::string::npos && c != '-' && c != ' '){
		  return false;	
		}
	}
	return true;
}

std::string 
Msa :: getCol(int col) const
{
  std::string column;
	for (int i(0); i < nseq; ++i){
		column.push_back(mali_seq[i][col]);
	}
	return column;
}


/**************************************************************
 * fitToAlphabet(string alph1) transforms symbols 
 * if a symbol from msa is not in alph1 then it is a gap
 **************************************************************/
void
Msa :: fitToAlphabet(const std::string & alph1){
	array<bool,256> allowed;
	allowed.fill(false);
	for (size_t i(0); i < alph1.size(); ++i){
		allowed[static_cast<unsigned char>(alph1[i])] = true;
	}

	std::string removed_symbols;
	array<bool,256> seen_removed;
	seen_removed.fill(false);
	for (int i(0); i < nseq; ++i){
		for (int j(0); j < ncol; ++j){
			const char symbol = mali_seq[i][j];
			if (symbol == '-' || symbol == ' '){
				continue;
			}
			const unsigned char c = static_cast<unsigned char>(symbol);
			if (!allowed[c]){
				mali_seq[i][j] = '-';
				if (!seen_removed[c]){
					seen_removed[c] = true;
					removed_symbols.push_back(symbol);
				}
			}
		}
	}

	for (size_t i(0); i < removed_symbols.size(); ++i){
		const char symbol = removed_symbols[i];
		const size_t pos = alphabet.find(symbol);
		if (pos != std::string::npos){
			alphabet.erase(alphabet.begin() + pos);
		}
		for (int col(0); col < ncol; ++col){
			const size_t aa_pos = aa_type_list[col].find(symbol);
			if (aa_pos != std::string::npos){
				aa_type_list[col].erase(aa_type_list[col].begin() + aa_pos);
				nb_type[col]--;
			}
		}
	}
	rebuildAlphaIndex();
}


/**************************************************************
 * printBasic() prints basic information in output
 *
 **************************************************************/
void
Msa :: printBasic(){
	std::string dictionary = "ARNDCQEGHILKMFPSTWYV-";
	std::vector<int> counts(dictionary.size(), 0);
	std::string out_name = Options::Get().output_fname;
	out_name = out_name.substr(0,out_name.find('.')) + ".aa_count";
	ofstream file(out_name.c_str());
	if (!file.is_open()){
		throw std::runtime_error("Cannot open file " + out_name);
	}
	for (int a(0); a < static_cast<int>(dictionary.size()); a++) {
		file << dictionary[a] << " ";
	}
	file << "\n";
	for (int col(0); col < ncol; col++){
		for (int seq(0); seq < nseq; seq++){
			int pos = static_cast<int>(dictionary.find(mali_seq[seq][col]));
			if (pos < static_cast<int>(dictionary.size())){
				counts[pos]++;
			} else {
				cerr << mali_seq[seq][col] << " is not in the dictionary\n";
			}
		}
		for (int a(0); a < static_cast<int>(dictionary.size()); a++) {
			file << counts[a] << " ";
			counts[a] = 0;
		}
		file << "\n";
	}
	file.close();
}


/**************************************************************
 * getSeqWeights() calculates the weight of each sequence
 * in the multiple alignment, following Henikoff & Henikoff (1994):
 *   w_i = \frac{1}{L}\sum_{x=1}^{L}\frac{1}{k_x n_{x,i}}
 * with
 *   L    = number of columns
 *   k_x  = number of distinct symbol types in column x (= nb_type[x])
 *   n_{x,i} = number of sequences sharing the same symbol as
 *             sequence i in column x
 *
 * Same formula as the calcSeqWeight() previously duplicated in
 * wentropy.cpp, trident.cpp and jensen.cpp, but computed for all
 * sequences at once in O(nseq*ncol) instead of O(nseq^2*ncol):
 * for each column, the occurrence count of every symbol (n_{x,a})
 * is precomputed once, instead of re-scanning the whole column
 * for every sequence. The result is cached: repeated calls (e.g.
 * from several statistics) cost nothing after the first one.
 **************************************************************/
const std::vector<float> &
Msa :: getSeqWeights(){
	if (seq_weight_computed){
		return seq_weight;
	}
	
	seq_weight = std::vector<float>(nseq, 0.0);
	std::vector<int> col_count(alphabet.size(), 0);
	
	for (int col(0); col < ncol; ++col){
		fill(col_count.begin(), col_count.end(), 0);
		for (int seq(0); seq < nseq; ++seq){
			col_count[getAaPos(mali_seq[seq][col])]++;
		}
		int k = nb_type[col];
		for (int seq(0); seq < nseq; ++seq){
			int n = col_count[getAaPos(mali_seq[seq][col])];
			seq_weight[seq] += 1.0 / (float) (n * k);
		}
	}
	for (int seq(0); seq < nseq; ++seq){
		seq_weight[seq] /= static_cast<float>(ncol);
	}
	
	seq_weight_computed = true;
	return seq_weight;
}

