/*
 *  wentropy.h
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 10/11/10.
 *
 */

#ifndef __WENTROPY_H__
#define __WENTROPY_H__

#include "statistic.h"

class WEntStat  : public Statistic 
{
private:
	int ncol;
	int nseq;
  float ** proba;						/**< probabilities of amino acids in each column (size = 20 * nb columns) */
	vector<float> seq_weight;	/**< Weight of each sequence in the msa (size = nb sequences) */
	vector<float> col_cons;		/**< Conservation score of each column (size = nb columns) */

	float calcSeqWeight(Msa & msa, int i);
	
public:
	~WEntStat();
	void calculateStatistic(Msa & msa);
};

#endif

