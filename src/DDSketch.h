/********************************************************/
/* AFQN Algorithm                                       */
/* Approximate Fast Qn in streaming                     */
/*                                                      */
/* Coded by Catiuscia Melle                             */
/*                                                      */
/* April 8, 2021                                        */
/*                                                      */
/* This code accompanies the paper                      */
/* AFQN: Approximate Qn Estimation in Data Streams      */
/*                                                      */
/* By: I. Epicoco, C. Melle, M. Cafaro and  M. Pulimeno */
/*                                                      */
/********************************************************/



#ifndef __DDSKETCH_H__
#define __DDSKETCH_H__

#include "Utility.h"
#include "IIS.h"

const int MIN_KEY = pow(2,30);                  


//****** ****** ****** ****** ****** ************ Utility functions


int getKeyFor(double value, double gamma, double logG);


double getCurrentAlpha(double alpha);


double getCurrentGamma(double alpha);


double getCurrentLogG(double gamma);


double getQuantileFraction(int kth, int I);


void debugSketch(std::map< int, int>& mySketch);



double estimateQ(std::map<int, int>& Sketch, double q, double gamma, int n);

void logQuantiles(FILE *fp, std::map< int,int>&mySketch, int collapses, double gamma, double *exactDiffs, int len);



//****** ****** ****** ****** ****** ************ Uniform Collapse

int performCollapse(std::map<int, int>& Sketch, int sketchBound, double *currentAlpha, double *currentGamma, double *currentLogG, int *SketchSize);

//****** ****** ****** ****** ****** ************ Sketch Filling (with s(s-1)/2 differences)

int fillSketch(int pos, double *window, double gamma, double LogG, std::map<int, int>& Sketch);


//****** ****** ****** ****** ****** ************ Sketch Updating


int updateSketch(double old_item, double new_item, double *Pwindow, int s, std::map<int,int>& Sketch, double gamma, double logGamma, int ndiffs);


int decreaseBinCount(int key, std::map<int, int>& sketch);


void updateSynopsis(double old_item, double new_item, double *Pwindow, int s, std::map<int,int>& Sketch, double gamma, double logGamma);


#endif //__DDSKETCH_H__

