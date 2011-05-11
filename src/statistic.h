/* Copyright (c) 2010 Guillaume Collet
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

#ifndef __STATISTIC_H__
#define __STATISTIC_H__

#include <vector>
#include <string>
#include <fstream>

#include "msa.h"
#include "options.h"
#include "factory.h"

using namespace std;

class Statistic 
{
public:
	Statistic(){};
	virtual ~Statistic(){};
	virtual void calculate(Msa & msa){};
	virtual void print(Msa & msa){};
};

class StatisticFactory : public Factory<Statistic>{};

void AddAllStatistics();

class Stat1D : public Statistic {
protected:
	vector<float> col_stat; /**< vector to store columns statistics */
	
public:
	virtual ~Stat1D(){};
	virtual void calculate(Msa & msa){};
	void print(Msa & msa){
		ofstream file(Options::Get().output_name.c_str());
		if (!file.is_open()){
			cerr << "Cannot open file " << Options::Get().output_name << "\n";
			exit(0);
		}
		if (Options::Get().global){
			float total = 0.0;
			for (int col(0); col < col_stat.size(); ++col){
				total += col_stat[col];
			}
			file << total / col_stat.size() << "\n"; 
		} else {
			for (int col(0); col < col_stat.size(); ++col){
				file << col + 1 << "\t" << col_stat[col] << "\n";
			}
		}
		file.close();
	};
};

class Stat2D : public Statistic {
protected:
	vector< vector<float> > cor_stat; /**< vector to store pairs of columns statistics */
	
public:
	virtual ~Stat2D(){};
	virtual void calculate(Msa & msa){};
	void print(Msa & msa){
		ofstream file(Options::Get().output_name.c_str());
		if (!file.is_open()){
			cerr << "Cannot open file " << Options::Get().output_name << "\n";
			exit(0);
		}
		for  (int x(0); x < cor_stat.size() - 1; ++x) {
			for (int y(0); y < cor_stat.size(); ++y) {
				if (y > x){
					file << cor_stat[x][y] << "\t";
				} else {
					file << 0.0 << "\t";
				}
			}
			file << "\n";
		}
		file.close();
	};
};

#endif

