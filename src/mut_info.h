/*
 *  mut_info.h
 *  MstatX++
 *
 *  Created by Guillaume Collet on 27/11/10.
 *
 */


#ifndef __MUT_INFO_H__
#define __MUT_INFO_H__

#include "statistic.h"

class MIStat  : public Statistic 
{
private:
	int ncol;				/**< number of columns in the multiple alignement */
	int npair;			/**< number of pairs of columns (npair = ncol * (ncol-1) / 2 */
	int nseq;				/**< number of sequences in the multiple alignement */
	int alph_size;	/**< size of the alphabet */
  float ** f_i;		/**< weighted frequencies of amino acids (size = alph_size * ncol) */
	float *** f_ij;	/**< weighted frequencies of pairs of amino acids (size = alph_size * alph_size * npair) */
	
	vector<float> seq_weight;		/**< Weight of each sequence in the msa (size = nb sequences) */
	vector<float> mutual_info;	/**< Mutual information in each pair of columns (size = npair) */
	
	float calcSeqWeight(Msa & msa, int i);
	
public:
	~MIStat();
	void calculateStatistic(Msa & msa);
};

#endif

