/*
 *  statistic.h
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 09/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __STATISTIC_H__
#define __STATISTIC_H__

#include <string>

#include "msa.h"
#include "factory.h"

using namespace std;

class Statistic 
{
public:
	Statistic(){};
	virtual ~Statistic(){};
		
	virtual void calculateStatistic(Msa & msa){};
};


class StatisticFactory : public Factory<Statistic>{};

void AddAllStatistics();


#endif