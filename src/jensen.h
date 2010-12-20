/*
 *  jensen.h
 *  MstatX
 *
 *  Created by COLLET Guillaume on 20/12/10.
 *
 */

#ifndef __JENSEN_H__
#define __JENSEN_H__

#include "statistic.h"

class JensenStat  : public Statistic 
{
private:
	int ncol;
	int nseq;
  double ** proba;						/**< probabilities of amino acids in each column (size = 20 * nb columns) */
	vector<double> seq_weight;	/**< Weight of each sequence in the msa (size = nb sequences) */
	vector<double> col_cons;		/**< Conservation score of each column (size = nb columns) */
	
	double calcSeqWeight(Msa & msa, int i);
	
public:
	~JensenStat();
	void calculateStatistic(Msa & msa);
};

#endif