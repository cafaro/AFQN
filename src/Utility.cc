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



#include "Utility.h"
#include <cstring>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern double NULLBOUND;
extern char VERSION[];  

// ******************************************************* DEBUG LOG

void logStartup(int s, int sketchBound, long N, int I, int kth, double quantile, int diff_frac, double currentAlpha, double currentGamma,  double QnScale) {
    
    std::cout << "\n\tApproximate Online Qn estimator, version "<< VERSION << std::endl;
    std::cout << "\tWindow size: " << s << std::endl; 
    std::cout << "\tLen of stream to process: " << N << std::endl;
    std::cout << "\tQn correction factor = " << QnScale << std::endl;
    std::cout << "\tInset size: " << I << std::endl; 
    std::cout << "\tk-th: "<< kth << std::endl;
    std::cout <<"\tq: " << quantile << std::endl; 
    std::cout << "\tSketch bound: " << sketchBound << std::endl;
    std::cout <<"\tFraction of differences computed at each iteration: " << diff_frac << std::endl; 
    std::cout << "\tInitial alpha = " << currentAlpha << ", gamma " << currentGamma << ", bound for near-0 values: " << NULLBOUND << std::endl; 
    std::cout<< "\n"<< std::endl;
    std::cout <<  "\t EXACT MEDIAN ESTIMATION\n";

    #ifdef TEST
        std::cout << "\t -----> TESTING MODE <-----\n\n";
    #elif DEBUG
        std::cout <<  "\t ***** DEBUG MODE *****\n\n";
    #elif VERIFY
        std::cout <<  "\t ***** VERIFY MODE *****\n\n";
    #else
        std::cout <<  "\t ++++++ CHECK MODE ++++++\n\n";
    #endif
}



void logWindow(double *window, int s, char type) {
    std::cout << "Window "<< type << "\t[ " ;
    for(int j = 0; j < s; ++j){
        std::cout << window[j] << " ";
    }
    std::cout << "] "<< std::endl;
}



void logDiff(std::vector<double>& diffs) {
    #ifdef DEBUG
    std::cout << "Differences\t[ " ;
    for (const auto& element : diffs) {
        std::cout << element<< " ";
    }
    std::cout << "] "<< std::endl;
    #endif
}

void debugDifferences(double diff, int key, double alpha, double gamma, double logG) {
    #ifdef DEBUG
    std::cout << "\t Computed difference " << diff << ", index: " << key;
    std::cout << ", [current (α, 𝛄, LogG) " << alpha;
    std::cout << ", " << gamma;
    std::cout << ", " << logG << " ]" << std::endl;
    #endif
}




// ******************************************************* EXACT COMPUTATION

double getExactKth(std::vector<double>& array, int kth) {
    std::nth_element(array.begin(), array.begin() + (kth-1), array.end());
    return array[kth-1];
}


double checkApproximationError(double estimate, double exact, std::string msg) {
    double error = std::abs((estimate-exact)/exact);  // approximation relative error  
    
    #ifdef DEBUG
        std::cout << msg << ": Exact value = " << exact;
        std::cout << " \tApproximate value = " << estimate;
        std::cout << " \tRelative Error (α) = " << error << std::endl;
    #endif

    return error;
}

// ******************************************************* Command Line Options

void initDefault(int *window_size, int *sketch_bound, int *streamLen, double *initial_alpha){
	(*window_size) = DEFAULT_WINDOW_SIZE;
	(*sketch_bound) = 2*(DEFAULT_WINDOW_SIZE);
	(*streamLen) = STREAMLEN;
	(*initial_alpha) = Alpha_0;
}


void printUsage(char *msg) {
    std::cerr << "Usage: " << msg << " {[-f path-to-file] | [-d distribution_type] [-x distribution_param] [-y distribution_param]} ";
    std::cerr << "[-s window_size] ";
    std::cerr << "[ -n max_stream_len ] [ -a initial_alpha ] [-b max_sketch_bound]\n\n" << std::endl;
    
    std::cerr << " -n is the len of the stream for the online phase (total items N = n+s)\n";
    std::cerr << " -d can be: \n";
    std::cerr << " : 1 Uniform distribution, with params [a:b] given by -x and -y options\n";
    std::cerr << " : 2 Exponential distribution, with params [λ] given by -x option\n";
    std::cerr << " : 3 Normal distribution, with params [µ:σ] given by -x and -y options\n";
    std::cerr << "\n";
}



void parseCommandLine(int argc, char *argv[], int *window_size, int *sketch_bound, double *initial_alpha, int *diff_fraction, Counters *stats){
	
    if (!stats) {
        fprintf(stderr, "ERROR: Counters variable is NULL\n");
        exit(1);
    }
    
    (*window_size) = DEFAULT_WINDOW_SIZE;
    (*sketch_bound) = 0; 
    (*initial_alpha) = Alpha_0;
    stats->MaxStreamLen = STREAMLEN;
    (*diff_fraction) = -1;

    bool file_flag = false;
    int t=-1;
    int c=0;
    while ( (c = getopt(argc, argv, "f:s:b:a:n:t:")) != -1) {
        
        switch (c) {
            
            case 'f':
                if (strlen(optarg) <= FSIZE) {
                    stats->filename = strndup(optarg, strlen(optarg));
                }
                file_flag = true;
                break;

            case 'n':
                stats->MaxStreamLen = strtol(optarg, NULL, 10);
                break;

            case 's':
				(*window_size) = atoi(optarg);
                break;        

            case 'b':
                (*sketch_bound) = atoi(optarg);
                break;
            
            case 'a':
                (*initial_alpha) = strtod(optarg, NULL); 
                break;

            case 't':
                t = atoi(optarg);
                break;

            default:
                fprintf(stderr, "?? getopt returned character code 0%o ??\n", c);
                break;
        }// switch
    }//wend

	#ifdef TEST
        if (!file_flag){
            std::cerr << "ERROR: a file with input data MUST be provided, -f path-to-file " << std::endl;
            exit(1);
        }
    #endif

    if ((*sketch_bound) == 0) {
        (*sketch_bound) = 2 * (*window_size);
    }

    if (t == -1) {
        t = 1; 
    } 
    (*diff_fraction) = std::ceil((double) ((*window_size) - 1 )/(t*1.0));

}



int checkCommandLineConfiguration(int argc, char *argv[], int *window_size, int *sketch_bound, double *initial_alpha, Counters *stats){
	
    int invalidRes = 1;

    if (!stats) {
        fprintf(stderr, "ERROR: Counters variable is NULL\n");
        return invalidRes;
    }
    
    (*sketch_bound) = 0; 
    (*window_size) = 0;         
    (*initial_alpha) = 0.0;     
    stats->streamLen = 0;
    stats->MaxStreamLen = 0;    

    bool file_flag = false;
    
    int distrtype = 0;
    double xparam = 0.0, yparam = 0.0;
    bool dist_flag = false;
    
    int c=0;
    while ( (c = getopt(argc, argv, "f:s:b:a:n:d:x:y:")) != -1) 
    {
        
        switch (c) 
        {
            
            case 'f':
                if (strlen(optarg) <= FSIZE) {
                    stats->filename = strndup(optarg, strlen(optarg));
                }
                file_flag = true;
                break;

            case 'n':
                stats->streamLen = strtol(optarg, NULL, 10);
                break;

            case 's':
				(*window_size) = atoi(optarg);
                break;        

            case 'b':
                (*sketch_bound) = atoi(optarg);
                break;
            
            case 'a':
                (*initial_alpha) = strtod(optarg, NULL); 
                break;

            case 'd':
                distrtype = atoi(optarg);
                dist_flag = true;
                break;
            
            case 'x':
                xparam = strtod(optarg, NULL); 
                break;

            case 'y':
                yparam = strtod(optarg, NULL); 
                break;

            default:
                fprintf(stderr, "?? getopt returned character code 0%o ??\n", c);
                break;
        }// switch
        
    }//wend getopt()

    if (*initial_alpha <= 0.0) {
        fprintf(stderr, "ERROR: α param not defined\n");
        return invalidRes;
    }

    if (! *window_size) {
        fprintf(stderr, "ERROR: window size not defined\n");
        return invalidRes;
    }

    if (! *sketch_bound) {
        fprintf(stderr, "ATTENTION: sketch bound not defined: setting on behalf of the window size\n");
        (*sketch_bound) = 2 * (*window_size);
    }

    if (!stats->streamLen){
        fprintf(stderr, "ERROR: total stream len N is equal to: s+n. You must provide -n\n");
        return invalidRes;
    }

    stats->MaxStreamLen = stats->streamLen + (*window_size);

    if (!file_flag && !dist_flag) {
        fprintf(stderr, "ERROR: at least an input file or a distribution type MUST be provided, -f or -d options\n");
        return invalidRes;
    }

    if (file_flag && dist_flag) {
        fprintf(stderr, "ERROR: you must provide or -f for an input file or -d for an input distribution, not both\n"); 
        return invalidRes;
    }


    #ifdef TEST
        
        if (!file_flag) {
            fprintf(stderr, "ERROR: in TEST Mode we want input from a file\n");;
            return invalidRes;
        }
    
    #else

        if (dist_flag) 
        {
            if (distrtype < 1 || distrtype > 3){
                fprintf(stderr, "ERROR: unrecognized distribution type (can be 1 or 2 or 3)\n");
                return invalidRes;
            }

            if (distrtype == 1 && ((xparam==0.0 && yparam==0.0) || (xparam >= yparam)) ) {
                fprintf(stderr, "ERROR: incorrect setting the range [a,b) for Uniform distribution\n");
                return invalidRes;
            }

            if (distrtype == 2 && xparam==0.0) {
                fprintf(stderr, "ERROR: incorrect setting the λ value for Exponential distribution\n");
                return invalidRes;
            }

            if (distrtype == 3 && (xparam==0.0 && yparam==0.0) ) {
                fprintf(stderr, "ERROR: incorrect setting mean (μ) and stddev (σ) for Normal distribution\n");
                return invalidRes;
            }

            stats->dtype = distrtype;
            stats->xparam = xparam;
            stats->yparam = yparam;
        }//fi dist_flag

    #endif

return 0;
}



void bufferStreamFromFile(Counters *stats) {

    if (stats && stats->filename) {

        FILE *fp = fopen(stats->filename, "r");
        if (fp == NULL) {
            fprintf(stderr,"Error opening %s\n", stats->filename);
            exit(1);
        }

        stats->item_points = (double *)malloc( sizeof(double) * stats->MaxStreamLen); 
        
        char *line = NULL;
        size_t dim = 0;
        long idx = 0;
        while( (getline(&line, &dim, fp)) != -1 && idx < stats->MaxStreamLen) {
            stats->item_points[idx] = strtod(line, NULL);
            ++idx;
        }//wend
        
        if (line)
            free(line);
        fclose(fp);

    }//fi
}



void initOutliersStats(Counters *stats) {

    stats->filename = NULL;
    stats->dtype = 0;
    stats->xparam = 0.0;
    stats->yparam = 0.0;

    stats->outlierFile = NULL;                              
    stats->inlierFile = NULL;    

    stats->approx_out_count = 0;
    stats->approx_in_count = 0;

    stats->streamLen = 0;
    stats->MaxStreamLen = 0;
    
    stats->item_points = NULL;

    stats->fpO = stats->fpI = NULL;
    
    #ifndef TEST
        stats->fpExactO = NULL;
        stats->fpExactI = NULL;

        stats->exact_out_count = 0;
        stats->exact_in_count = 0;
        
        stats->exac_outF = NULL;
        stats->exac_inF = NULL;
    #else
        stats->outliersBuffer = NULL;
        stats->inliersBuffer = NULL; 
    #endif
}

void initResultFilename(Counters *stats, int window_size, int tsize) {
  
    std::string name = stats->filename;
    std::size_t pos = name.find_last_of("/"); 
    std::string distr = name.substr(pos+1, 4); 
    
    std::string n1 = "";
    std::string n2 = "";

    n1 += "./" + distr + "-Outlier-" + std::to_string(window_size) + "-" + std::to_string(tsize) + ".csv";
    n2 += "./" + distr + "-Inlier-" + std::to_string(window_size) + "-" + std::to_string(tsize) + ".csv";
    
    stats->outlierFile = new char[n1.length()+1];
    std::strcpy(stats->outlierFile, n1.c_str());

    stats->inlierFile = new char[n2.length()+1];
    std::strcpy(stats->inlierFile, n2.c_str());
}


void initExactFilename(Counters *stats, int window_size, int tsize) {

    #ifndef TEST
        if (stats!= NULL && stats->filename != NULL) {
            std::string name = stats->filename;
            std::size_t pos = name.find_last_of("/"); 
            std::string distr = name.substr(pos+1, 4); 
    
            std::string n1 = "";
            std::string n2 = "";

            n1 += "./" + distr + "-ExactOutlier-" + std::to_string(window_size) + "-" + std::to_string(tsize) + ".csv";
            n2 += "./" + distr + "-ExactInlier-" + std::to_string(window_size) + "-" + std::to_string(tsize) + ".csv";

            stats->exac_outF = new char[n1.length()+1];
            std::strcpy(stats->exac_outF, n1.c_str());
                
            stats->exac_inF = new char[n2.length()+1];
            std::strcpy(stats->exac_inF, n2.c_str());
        }//fi
    #endif  
}


void openLog(Counters *stats) {

    if (stats->outlierFile) {
        stats->fpO = fopen(stats->outlierFile, "w");
        if (stats->fpO == NULL) {
            fprintf(stderr, "Error opening %s\n", stats->outlierFile);
            exit(1);
        } 
        fprintf(stats->fpO, "%s,%s,%s,%s,%s,%s,,%s\n", "seqNo", "item", "Median", "Q1", "z-score", "collapse", "alpha");
    }

    if (stats->inlierFile){
        stats->fpI = fopen(stats->inlierFile, "w");
        if (stats->fpI == NULL) {
            fprintf(stderr, "Error opening %s\n", stats->inlierFile);
            exit(1);
        } 
        fprintf(stats->fpI, "%s,%s,%s,%s,%s,%s,,%s\n", "seqNo", "item", "Median", "Q1", "z-score", "collapse", "alpha");
    }

    #ifdef TEST
        stats->outliersBuffer = (Item *)malloc(sizeof(Item)*stats->MaxStreamLen);                           
        stats->inliersBuffer  = (Item *)malloc(sizeof(Item)*stats->MaxStreamLen); 
    #else 
        stats->fpExactO = NULL;
        stats->fpExactI = NULL;            

        if (stats->exac_outF){
            stats->fpExactO = fopen(stats->exac_outF, "w");
            if (stats->fpExactO == NULL) {
                fprintf(stderr, "Error opening %s\n",stats->exac_outF);
                exit(1);
            } 
            fprintf(stats->fpExactO, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "seqNo", "item", "Median", "K-th", "Q1", "relErr", "Qn", "z-score", "collapse", "#bins", "alpha");
        }

        if (stats->exac_inF){
            stats->fpExactI = fopen(stats->exac_inF, "w");
            if (stats->fpExactI == NULL) {
                fprintf(stderr, "Error opening %s\n",stats->exac_inF);
                exit(1);
            } 
            fprintf(stats->fpExactI, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "seqNo", "item", "Median", "K-th", "Q1", "relErr", "Qn", "z-score", "collapse", "#bins", "alpha");
        }
    #endif 
}


void closeLog(Counters *stats) {

    #ifdef TEST
        if (stats->outliersBuffer){
            for(long u = 0; u<(stats->approx_out_count); ++u) {
                fprintf(stats->fpO, "%ld,%.6f\n", stats->outliersBuffer[u].seq, stats->outliersBuffer[u].middle);
            } 
            free(stats->outliersBuffer);
        }
            
        if (stats->inliersBuffer){
            for(long u = 0; u<(stats->approx_in_count); ++u) {
                fprintf(stats->fpI, "%ld,%.6f\n", stats->inliersBuffer[u].seq, stats->inliersBuffer[u].middle);
            } 
            free(stats->inliersBuffer);
        }
    #else
        if (stats->fpExactO) {
            fclose(stats->fpExactO);
        }
        
        if (stats->fpExactI) {
            fclose(stats->fpExactI);
        }
    #endif

    if (stats->fpO) {
        fclose(stats->fpO);
    }
    
    if (stats->fpI) {
        fclose(stats->fpI);
    }
}


void destroyOutliersStats(Counters *stats) {
   
    if (stats) {
    
        if (stats->filename) {
            free(stats->filename);
        }

        //filenames are allocated with new
        if (stats->outlierFile){
            delete stats->outlierFile;
        }

        if (stats->inlierFile){
            delete stats->inlierFile;
        }

        if (stats->item_points){
            free(stats->item_points);
        }
 
        #ifndef TEST
            if (stats->exac_outF) {
                delete stats->exac_outF;
            }
            
            if (stats->exac_inF) {
                delete stats->exac_inF;
            }
        #endif
    }//fi
}

// ******************************************************* OUTLIERNESS

double getQnScaleFactor(int n, double scalingFactor) {

	double dn;
	if (n <= 9) {
		if (n == 2) {
			dn = .399;
		} else if (n == 3) {
			dn = .994;
		} else if (n == 4) {
			dn = .512;
		} else if (n == 5) {
			dn = .844;
		} else if (n == 6) {
			dn = .611;
		} else if (n == 7) {
			dn = .857;
		} else if (n == 8) {
			dn = .669;
		} else {
			dn = .872;
		}
	} else {
		if (n % 2 == 1) {
			dn = n / (n + 1.4);
		} else {
			dn = n / (n + 3.8);
		}
	}//fi (n<9)
    return  (dn * scalingFactor);
}



void setQnValue(Counters *stats, int window_size) {

    // as per Rousseeuw and Croux (paper 1992): 
    //QFactor = 2.2219;                  
    stats->QnScale = getQnScaleFactor(window_size, QFactor);          //Qn correction factor
}



void checkForOutlier(double middle, long seqNo, double median, double Q1, Counters *stats, double alpha, int collapse, int bins) {
    
    double t = 3.0;
    double zscore = (fabs(middle - median) - (t * stats->QnScale * Q1));
    
    if ( zscore > 0) {
        
		#ifndef TEST
			fprintf(stats->fpO, "%ld,%.6f,%.6f,%.6f,%.6f,%d,%d,%.6f\n", seqNo, middle, median, Q1, zscore, collapse, bins, alpha);
        #else
            stats->outliersBuffer[stats->approx_out_count].seq = seqNo;
            stats->outliersBuffer[stats->approx_out_count].middle = middle;    
		#endif
        
        ++stats->approx_out_count;
	} else {
        
        #ifndef TEST
			fprintf(stats->fpI, "%ld,%.6f,%.6f,%.6f,%.6f,%d,%d,%.6f\n", seqNo, middle, median, Q1, zscore, collapse, bins, alpha);
        #else
            stats->inliersBuffer[stats->approx_in_count].seq = seqNo;
            stats->inliersBuffer[stats->approx_in_count].middle = middle;     
		#endif
        
        ++stats->approx_in_count;
    }//fi check test
}




void exactOutlier(double middle, long seqNo, double exactM, double exactK, Counters *stats, double apprK, double errQ, double alpha, int collapse, int bins) {

    double t = 3.0;
    double Qn = stats->QnScale * exactK;
    
    double zscore = (fabs(middle - exactM) - (t*Qn));

    if (zscore > 0) {
            
        #ifndef TEST
            ++stats->exact_out_count;

            fprintf(stats->fpExactO, "%ld,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%.6f\n", seqNo, middle, exactM, 
            exactK, apprK, errQ, Qn, zscore, collapse, bins, alpha);
        #endif
	} else {

        #ifndef TEST
            ++stats->exact_in_count;
		
        	fprintf(stats->fpExactI, "%ld,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%.6f\n", seqNo, middle, exactM, 
            exactK, apprK, errQ, Qn, zscore, collapse, bins, alpha);
        #endif
    }//fi check test
}




void OutlierTest(double middle, long seqNo, double median, double Q1, Counters *stats, int collapse, double alpha) {
    
    double t = 3.0;
	double zscore = (fabs(middle - median) - (t*(stats->QnScale)*Q1));

    if (zscore > 0) {

    	#ifdef TEST
            stats->outliersBuffer[stats->approx_out_count].seq = seqNo;
            stats->outliersBuffer[stats->approx_out_count].middle = middle;    
        #else
    		fprintf(stats->fpO, "%ld,%.6f,%.6f,%.6f,%.6f,%d,%.6f\n", seqNo , middle, median, Q1, zscore, collapse, alpha);
		#endif
        ++stats->approx_out_count;
	} else {

        #ifdef TEST
            stats->inliersBuffer[stats->approx_in_count].seq = seqNo;
            stats->inliersBuffer[stats->approx_in_count].middle = middle;    
        #else 
            fprintf(stats->fpI, "%ld,%.6f,%.6f,%.6f,%.6f,%d,%.6f\n", seqNo , middle, median, Q1, zscore, collapse, alpha);
        #endif
        ++stats->approx_in_count;
    }//fi check test
}



void startTimer(Timer *t) {
    if (t) {
        gettimeofday(&t->start, NULL);
    }
}

void stopTimer(Timer *t) {
    if (t) {
        gettimeofday(&t->end, NULL);
    }
}

double getElapsedMilliSecs(Timer *t) {
    
    if (t) {
        double seconds = (t->end.tv_sec - t->start.tv_sec);
        double microsec = (t->end.tv_usec - t->start.tv_usec);
        return (1000*seconds+((double)microsec/1000.0));
    } else {
        return -1;
    }
}


double getElapsedSeconds(Timer *t){
    
    if (t) {
        double seconds = (t->end.tv_sec - t->start.tv_sec);
        double microsec = (t->end.tv_usec - t->start.tv_usec);
        return (seconds+((double)microsec/1000000.0));
    } else {
        return -1;
    }
}

