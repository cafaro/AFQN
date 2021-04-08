/*
 * (Some sample C code for the quickselect algorithm, taken from Numerical Recipes in C):
 * http://www.stat.cmu.edu/~ryantibs/median/quickselect.c
 * 
 */

#ifndef __QSELECT_H__
#define __QSELECT_H__

#include <vector>

#include <cmath>
#include <iterator>
#include <utility>

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;




double quickselect(double *data, long len, long pos);

double selectQuickly(std::vector<double> &data, long len, long pos);

#endif //__QSELECT_H__
