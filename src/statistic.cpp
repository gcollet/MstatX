/*
 *  statistic.cpp
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 09/11/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "statistic.h"
#include "wentropy.h"
#include "mut_info.h"
#include "trident.h"
#include "kabat.h"
#include "mlc.h"

void AddAllStatistics()
{
	StatisticFactory::Add<MlcStat>("mlc");
	StatisticFactory::Add<TridStat>("trident");
	StatisticFactory::Add<WEntStat>("wentropy");
	StatisticFactory::Add<MIStat>("mutual");
	StatisticFactory::Add<KabatStat>("kabat");
}
