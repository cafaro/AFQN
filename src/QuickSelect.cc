/*
 * (Some sample C code for the quickselect algorithm, taken from Numerical Recipes in C):
 * http://www.stat.cmu.edu/~ryantibs/median/quickselect.c
 *
 */



#include "QuickSelect.h"


double quickselect(double *data, long len, long pos) {
  
  long i, ir, j, l, mid;
  double a, temp;

  l=0;
  ir=len-1;

  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && data[ir] < data[l]) {
	    SWAP(data[l],data[ir]);
      }
    return data[pos];
    }
    else 
    {
      mid=(l+ir) >> 1; 
      SWAP(data[mid],data[l+1]);
      if (data[l] > data[ir]) {
	    SWAP(data[l],data[ir]);
      }
      if (data[l+1] > data[ir]) {
	    SWAP(data[l+1],data[ir]);
      }
      if (data[l] > data[l+1]) {
	    SWAP(data[l],data[l+1]);
      }
      i=l+1; 
      j=ir;
      a=data[l+1];

      for (;;) { 
	    do i++; while (data[i] < a); 
	    do j--; while (data[j] > a); 
	    if (j < i) break; 
	    SWAP(data[i],data[j]);
      } 
      data[l+1]=data[j]; 
      data[j]=a;
      if (j >= pos) ir=j-1; 
      if (j <= pos) l=i;
    }
  }

}


double selectQuickly(std::vector<double> &data, long len, long pos) {
  
  long i, ir, j, l, mid;
  double a, temp;

  l=0;
  ir=len-1;

  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && data[ir] < data[l]) {
	    SWAP(data[l],data[ir]);
      }
    return data[pos];
    }
    else 
    {
      mid=(l+ir) >> 1; 
      SWAP(data[mid],data[l+1]);
      if (data[l] > data[ir]) {
	    SWAP(data[l],data[ir]);
      }
      if (data[l+1] > data[ir]) {
	    SWAP(data[l+1],data[ir]);
      }
      if (data[l] > data[l+1]) {
	    SWAP(data[l],data[l+1]);
      }
      i=l+1; 
      j=ir;
      a=data[l+1];

      for (;;) { 
	    do i++; while (data[i] < a); 
	    do j--; while (data[j] > a); 
	    if (j < i) break; 
	    SWAP(data[i],data[j]);
      } 
      data[l+1]=data[j]; 
      data[j]=a;
      if (j >= pos) ir=j-1; 
      if (j <= pos) l=i;
    }
  }

}
