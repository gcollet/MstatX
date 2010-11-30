/*
 *  trident.h
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 09/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __TRIDENT_H__
#define __TRIDENT_H__

#include "statistic.h"

class TridStat : public Statistic {
private:
	int ncol;
	int nseq;
	vector<float> seq_weight;	/**< Weight of each sequence in the msa (size = nb sequences) */
	vector<float> t;					/**< t(x) = Shannon entropy score  + Weighted sequence Score */
	vector<float> r;					/**< r(x) = Stereochemical score */
	vector<float> g;					/**< g(x) = Gap Score */
	vector<float> col_cons;		/**< Conservation score = (1-t(x))^a * (1- r(x))^b * (1-g(x))^c */
	
	float calcSeqWeight(Msa & msa, int i);
	float normVect(vector<float> vect);
	
public:
	void calculateStatistic(Msa & msa);
};

#endif