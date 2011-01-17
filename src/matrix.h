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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <string>
#include <iostream>

using namespace std;

template<class T>
class TriangularMatrix
{
private:
	int size;
	int width;
	T * data;
public:
	/* Constructor */
	TriangularMatrix():size(0),width(0){}
	TriangularMatrix(int width):width(width){
		size = width * (width + 1) / 2;
		data = new T[size];
	}
	
	/* Destructor */
	~TriangularMatrix(){
		delete [] data;
	}
	
	/* Initialization */
	void initMatrix(int w){
		width = w;
		size  = width * (width + 1) / 2;
		if (data != NULL){
			cerr << "data allocated in matrix are erased\n";
			free(data);
		}
		data = new T[size];
	}
	
	/* Setter */
	void set(int row, int col, T val){
	  int i = row;
		int j = col;
		if (i > j){
			i = col;
			j = row;
		}
		data[j * (j + 1) / 2 + i] = val;
	}
	
	/* Getter */
	T get(int row, int col){
	  int i = row;
		int j = col;
		if (i > j){
			i = col;
			j = row;
		}
		return data[j * (j + 1) / 2 + i];
	}
	
	/* Index */
	int index (int row, int col){
		int i = row;
		int j = col;
		if (i > j){
			i = col;
			j = row;
		}
		return j * (j + 1) / 2 + i;
	}
};

#endif
