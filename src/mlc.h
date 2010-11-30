/*
 *  mlc.h
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 09/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MLC_H__
#define __MLC_H__

#include <vector>

#include "statistic.h"

class MlcStat : public Statistic 
{
private:
	int ncol;
	float ** pair_mat;
  vector<float> mean;
	vector<float> stdev;
	
	float ** cor_mat;

public:
	~MlcStat();
	void calculateStatistic(Msa & msa);
};

#endif

