/*
 *  matrix.h
 *  MstatX
 *
 *  Created by COLLET Guillaume on 08/11/10.
 *
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
	TriangularMatrix():size(0),width(0),data(NULL){}
	TriangularMatrix(int width):width(width){
		size = width * (width + 1) / 2;
		
	}
	
	/* Destructor */
	~TriangularMatrix(){
		if ( data != NULL) {
		  free(data);	
		}
	}
	
	/* Initialization */
	void initMatrix(int w){
		width = w;
		size  = width * (width + 1) / 2;
		if (data != NULL){
			cerr << "data allocated in matrix are erased\n";
			free(data);
		}
		data = (T *) malloc(size * sizeof(T));
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
		return data[j * (j + 1) / 2 + i];
	}
};

#endif