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



#include "DDSketch.h"
#include "IIS.h"
#include "QuickSelect.h"
#include "Utility.h"

char VERSION[] = "AFQNv1";      
double NULLBOUND;               


int main(int argc, char *argv[]) {
        
    if (argc == 1) {
        printUsage(argv[0]);
        exit(1);
    }//fi argc
    
    // *********************** STREAM STATISTICS 
    
    Counters stats; 
    initOutliersStats(&stats);                      

    
    // *********************** PROCESSING CONSTRAINTS
    
    int s;                                         
    long N;                                        
    int sketchBound;                               
    double alpha;                                 
    int diff_fraction = 1;

    
    
    int isNotValid = checkCommandLineConfiguration(argc, argv, &s, &sketchBound, &alpha, &stats);   
    if (isNotValid) {
        std::cerr << " Command line configuration and options are not valid\n";
        return isNotValid;
    }

    // *********************** TIME (SLIDING) WINDOW
    
    double window[s];                              
    long seqNo[s];                                 
    long sLen = 0;                                 
    int pos = -1;                                  
    int middle_index = s/2;                        
    
    
    // *********************** (SLIDING) MEDIAN OF THE TIME WINDOW    
    int median_index = s/2;                          
    double Pwindow[s];                             
    int IIS_pos = -1;                              
    
    // *********************** Qn OF THE TIME WINDOW

    int h = s/2 + 1;                               
    int kth = h*(h-1)/2;                           
    
    setQnValue(&stats, s);                         
    int I = s*(s-1)/2;                             
    
    double quantile = getQuantileFraction(kth, I); 
    
    // *********************** INPUT STREAM and LOGS 
    bool filemode = true;                          
    
    #ifdef TEST  
        bufferStreamFromFile(&stats);               
        initResultFilename(&stats, s, sketchBound); 
    #else    

        // define the stream distribution params 
        std::uniform_real_distribution<double> udistribution(stats.xparam, stats.yparam);
        std::exponential_distribution<double> edistribution(stats.xparam);
        std::normal_distribution<double> ndistribution(stats.xparam, stats.yparam);
        
        std::default_random_engine generator;
        std::function<double()> randomizer;
        generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

        if (stats.filename != NULL) { 
    
            initResultFilename(&stats, s, sketchBound);
            initExactFilename(&stats, s, sketchBound);

        } else {
            filemode = false;  
            
            if (stats.dtype == 1) {
                stats.filename = strndup("Uniform", strlen("Uniform"));
                randomizer = std::bind(udistribution, generator);
            } else if (stats.dtype == 2) {
                stats.filename = strndup("Exponential", strlen("Exponential"));
                randomizer = std::bind(edistribution, generator);
            } else {
                stats.filename = strndup("Normal", strlen("Normal"));
                randomizer = std::bind(ndistribution, generator);
            }//fi dtype
            
            initResultFilename(&stats, s, sketchBound);
            initExactFilename(&stats, s, sketchBound);
        }//fi
    #endif

    openLog(&stats);

    #ifdef CMP
    Item *loggedPoints = (Item *)malloc(sizeof(Item)* stats.streamLen); 
    long pIdx = 0;
    #endif

    // *********************** SKETCH vars    
    std::map< int, int > Sketch;                          

    double currentAlpha = alpha;                          
    double currentGamma = getCurrentGamma(currentAlpha);  
    double currentLogG = getCurrentLogG(currentGamma);    
    
    NULLBOUND = pow(currentGamma, -MIN_KEY);              
    int Sketch_population = 0;                            
    int Sketch_size = 0;                                  
    int TotalCollapse = 0;                                

    double item;                                          
    double oldest_item;                                   
    double estimatedQ, exact_kth;                         
    double exact_M;                                       

    #ifndef TEST
        double new_diff, old_diff;                        
        double errQ;                                      
    
        std::vector<double> ExactDiffs(I,0);              
        int dd = 0;                                       

        char qfilename[FSIZE];
        std::string name = stats.filename;
        std::size_t from_p = name.find_last_of("/"); 
        std::string distr = name.substr(from_p+1, 4); 
        snprintf(qfilename, FSIZE-1, "./Quantiles-%s-%d-%d.csv", distr.c_str(), s, diff_fraction);
        
        FILE *qfile = fopen(qfilename, "w");
        if (qfile == NULL){
            fprintf(stderr,"Error opening %s\n", qfilename);
            exit(1);       
        }//fi
        fprintf(qfile, "Population,Bins,Collapses,EMin,Amin,err,index,EQ1,AQ1,err,index,EQ2,AQ2,err,index,EQ3,AQ3,err,index,EMax,AMax,err,index\n");
    #endif


    // ************************************ Starting processing
    logStartup(s, sketchBound, stats.MaxStreamLen, I, kth, quantile, diff_fraction, currentAlpha, currentGamma, stats.QnScale);   
    
    #ifdef TEST
        item = stats.item_points[sLen];
    #else
        if (!filemode) {
            item = randomizer();
        } else {
            item = stats.item_points[sLen];
        }
    #endif

    ++sLen;                     
    ++pos;                      
    window[pos] = item;         
    seqNo[pos] = sLen;                                            
    Pwindow[0] = item;          

    while (sLen < s) {
        #ifdef TEST
            item = stats.item_points[sLen];
        #else
            if (!filemode) {
                item = randomizer();
            } else {
                item = stats.item_points[sLen];
            }    
        #endif
        
        ++sLen;                                   
        ++pos;
        window[pos] = item;                        
        seqNo[pos] = sLen;
        
        IIS_pos = isort_v5(Pwindow, pos, item);    
        Sketch_population += fillSketch(pos, window, currentGamma, currentLogG, Sketch);
        TotalCollapse += performCollapse(Sketch, sketchBound, &currentAlpha, &currentGamma, &currentLogG, &Sketch_size);    
        
        #ifndef TEST
            for(int j = pos-1; j>=0; --j) { 
                new_diff = std::abs(window[pos]-window[j]);
                ExactDiffs[dd] = new_diff;          
                ++dd;
            }//for
        #endif

    }//wend
    
    exact_M = Pwindow[median_index];                                                              

    #ifndef TEST
        logQuantiles(qfile, Sketch, TotalCollapse, currentGamma, ExactDiffs.data(), I);
    #endif
    
    
    Timer onlineTime;        
    
    Timer medianT, quantileT, loopT;
    double medianS=0, quantileS=0, loopS=0, sketchS = 0;
    
    long countchecks = 0;                         
    
    startTimer(&onlineTime);
    for (long i = 0; i < stats.streamLen; i++) {
    
        #ifdef TEST
            item = stats.item_points[sLen];
        #else
            if (!filemode) {
                item = randomizer();
            } else {
                item = stats.item_points[sLen];
            }
        #endif
        
        ++sLen;
        pos = (pos+1)%s;                   
        oldest_item = window[pos];          
        window[pos] = item;                 
        seqNo[pos] = sLen;
        
        if (oldest_item != item)
        {
            updateSynopsis(oldest_item, item, Pwindow, s, Sketch, currentGamma, currentLogG);
            TotalCollapse += performCollapse(Sketch, sketchBound, &currentAlpha, &currentGamma, &currentLogG, &Sketch_size);
        
        }//fi

        #ifndef TEST
            for(int l=1; l<s; ++l) {                                                      
                int id = (pos+l)%s; 
                new_diff = std::abs(window[pos] - window[id]);
                old_diff = std::abs(oldest_item - window[id]);
                std::vector<double>::iterator it = std::find(ExactDiffs.begin(), ExactDiffs.end(), old_diff);
                if (it != ExactDiffs.end()){
                    *(it) = new_diff;
                }//fi
            }//for l
            logQuantiles(qfile, Sketch, TotalCollapse, currentGamma, ExactDiffs.data(), I);
        #endif
        
        exact_M = Pwindow[median_index];

        estimatedQ = estimateQ(Sketch, quantile, currentGamma, I);
        ++middle_index;                  

        #ifdef CMP        
            
            loggedPoints[pIdx].seq = seqNo[middle_index%s];
            loggedPoints[pIdx].middle = window[middle_index%s];
            loggedPoints[pIdx].median = exact_M;
            loggedPoints[pIdx].Qn = stats.QnScale * estimatedQ;
            loggedPoints[pIdx].collapses = TotalCollapse;
            loggedPoints[pIdx].alpha = currentAlpha;
            loggedPoints[pIdx].bins = Sketch.size();
            
            if ( (fabs(window[middle_index%s] - exact_M) - (3 * loggedPoints[pIdx].Qn)) > 0 ){
                loggedPoints[pIdx].isOutlier = 1;
                ++(stats.approx_out_count);
            }else{
                loggedPoints[pIdx].isOutlier = 0;
                ++(stats.approx_in_count);
            }//fi check

            ++pIdx;   
        #else
                                                           
            #ifdef TEST
                OutlierTest(window[middle_index%s], seqNo[middle_index%s], exact_M, estimatedQ, &stats, TotalCollapse, currentAlpha);
            #else
                exact_kth = quickselect(ExactDiffs.data(), I, kth-1 );
                
                errQ = checkApproximationError(estimatedQ, exact_kth, "k-th order stat");   
                
                checkForOutlier(window[middle_index%s], seqNo[middle_index%s], exact_M, estimatedQ, &stats, currentAlpha, TotalCollapse, Sketch_size);
            
                exactOutlier(window[middle_index%s], seqNo[middle_index%s], exact_M, exact_kth, &stats, estimatedQ, errQ, currentAlpha, TotalCollapse, Sketch_size);
            #endif
       #endif
        
        ++countchecks;                 
    }//for distribution len
    stopTimer(&onlineTime);

    double running_secs = (getElapsedMilliSecs(&onlineTime)/1000.0); 
    double update_per_sec = ((double)countchecks/running_secs);


    std::cout << "Processing "<< stats.filename << " ended" << std::endl;
    #ifndef TEST
        std::cout << "\nFound Exact Outliers "<< stats.exact_out_count << "  and " << stats.exact_in_count << " Exact Inliers in stream of length " << sLen << " items\n" << std::endl;
        std::cout << "\nFound Approximated Outliers "<< stats.approx_out_count << "  and " << stats.approx_in_count << " Approximated Inliers over StreamTotalLength " << N << std::endl;
        std::cout << "\nProcessing time (online phase only): "<< getElapsedMilliSecs(&onlineTime) << " ms " << std::endl;
        std::cout << "Collapse executed " << TotalCollapse << ", Final Alpha " << currentAlpha << ", Final Gamma " << currentGamma;
        std::cout << ", Final Bins " << Sketch.size() << std::endl;
    #else
        std::cerr << stats.filename << "," << countchecks << "," << s/2 << "," << running_secs << "," << update_per_sec;
        std::cerr << "," <<  stats.approx_out_count << "," << stats.approx_in_count;
        std::cerr << "," << alpha << "," << sketchBound;
        std::cerr << "," << TotalCollapse << "," << currentAlpha << "," << Sketch.size() << std::endl;
    #endif


    #ifdef CMP
        char *sub = strrchr(stats.filename, '/');
        char fname[FSIZE];
        
        int len = strlen(&sub[1]);
        char stripped[len-4];
        strncpy(stripped, &sub[1], len-4);
        stripped[len-4]='\0';
        
        snprintf(fname, FSIZE-1, "Results/%s-%d-%d.csv", stripped, s, sketchBound);
        
        FILE *logF = fopen(fname, "w");
        if (logF != NULL) {
        
            for(long u = 0; u<pIdx; ++u) {
                fprintf(logF, "%ld,%.6f,%.6f,%.6f,%d,%d,%.6f,%d\n", loggedPoints[u].seq, loggedPoints[u].middle, loggedPoints[u].median, loggedPoints[u].Qn, loggedPoints[u].isOutlier, loggedPoints[u].collapses,loggedPoints[u].alpha,loggedPoints[u].bins);
            }//for 

            fclose(logF);
            free(loggedPoints);
        }//fi
    #endif

    closeLog(&stats);
    destroyOutliersStats(&stats);

    #ifndef TEST
        fclose(qfile);
    #endif

    std::cout << "Processing ended!\n\n";
    return 0;
}

