/*
 *  kabat.h
 *  MstatX++
 *
 *  Created by COLLET Guillaume on 16/12/10.
 *
 */

#ifndef __KABAT_H__
#define __KABAT_H__

#include "statistic.h"

class KabatStat : public Statistic
{
private:
	vector<float> col_cons;		/**< Conservation score of each column (size = nb columns) */
	
public:
	void calculateStatistic(Msa & msa);
};

#endif