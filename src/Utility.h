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

#ifndef __UTILITY_H__
#define __UTILITY_H__


#include <iostream>  
#include <iomanip> 
#include <algorithm>
#include <random>
#include <cmath>
#include <iterator>
#include <map>
#include <vector>
#include <utility>

#include <stdio.h>
#include <sys/time.h>


const int DEFAULT_WINDOW_SIZE = 1001;           
const int STREAMLEN = 1001;     
const int PARSING_ERROR = 7;                    
const int FSIZE = 256;                          
const double QFactor = 2.2219;                  
const double Alpha_0 = 0.001;                   



typedef struct Item {
    long seq;           
    double middle;      
    
    double median;      
    double Qn;          

    int isOutlier;      

    int collapses;      
    double alpha;       
    int bins;           
} Item;



typedef struct LogPoint{
    int seqNo;
    double value;
    int isOutlier;
    double median;
    double qn;
    double zscore;
}LogPoint;



typedef struct PointDiff {
    int pos;
    
    double d1;
    double w1;
    
    double d2;
    double w2;

    double d3;
    double w3;
} PointDiff;


typedef struct Counters {
                            
    char *filename;             
    
    int dtype;                  
    double xparam;              
    double yparam;              

    char *outlierFile;          
    char *inlierFile;           

    FILE *fpO;                  
    FILE *fpI;                  
    
    double *item_points;        
    long streamLen;             
    long MaxStreamLen;          
    
    double QnScale;             
    
    int approx_out_count;       
    int approx_in_count;        

    #ifndef TEST
        FILE *fpExactO;         
        FILE *fpExactI;         
        char *exac_outF;        
        char *exac_inF;         
        int exact_out_count;    
        int exact_in_count;     
    #else  
        Item *outliersBuffer;   
        Item *inliersBuffer;    
    #endif 

} Counters;




// ******************** DEBUG VIEWS

void logStartup(int s, int sketchBound, long N, int I, int kth, double quantile, int diff_frac, double currentAlpha, double currentGamma,  double QnScale);



void logWindow(double *window, int s, char type);



void logDiff(std::vector<double>& diffs);



void debugDifferences(double diff, int key, double alpha, double gamma, double logG);


// ******************** EXACT QUANTILES

double getExactKth(std::vector<double>& array, int kth);


double checkApproximationError(double estimate, double exact, std::string msg);



// ******************** Management of command line arguments


void printUsage(char *msg);

void parseCommandLine(int argc, char *argv[], int *window_size, int *sketch_bound, double *initial_alpha, int *diff_fraction, Counters *stats);


int checkCommandLineConfiguration(int argc, char *argv[], int *window_size, int *sketch_bound, double *initial_alpha, Counters *stats);

// ******************** Files management 

void bufferStreamFromFile(Counters *stats);

void initOutliersStats(Counters *stats);

void openLog(Counters *stats);

void closeLog(Counters *stats);

void destroyOutliersStats(Counters *stats);

void initResultFilename(Counters *stats, int window_size, int tsize);

void initExactFilename(Counters *stats, int window_size, int tsize);

// ******************** Outlierness

double getQnScaleFactor(int n, double scalingFactor);
void setQnValue(Counters *stats, int window_size);

void OutlierTest(double middle, long seqNo, double median, double Q1, Counters *stats, int collapse, double alpha);

void exactOutlier(double middle, long seqNo, double exactM, double exactK, Counters *stats, double apprK, double errQ, double alpha, int collapse, int bins);
void checkForOutlier(double middle, long seqNo, double median, double Q1, Counters *stats, double alpha, int collapse, int bins);





// ****************** Time evaluation

typedef struct Timer{
    struct timeval start;
    struct timeval end;
}Timer;


void startTimer(Timer *t);
void stopTimer(Timer *t);

double getElapsedMilliSecs(Timer *t);
double getElapsedSeconds(Timer *t);


#endif //__UTILITY_H__
