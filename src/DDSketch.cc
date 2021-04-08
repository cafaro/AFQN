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
#include "QuickSelect.h"

extern double NULLBOUND;   

//****** ****** ****** ****** ****** ************ ************ ************ ************ ****** Utility functions

int getKeyFor(double value, double gamma, double logG) {
    
    if (value <= NULLBOUND) {
        #ifdef DEBUG
            std::cout << "\tgetKeyFor() "<< value << "\tDifference in Bucket 0 " << std::endl;
        #endif
        return -MIN_KEY; //bucket for 0
    }
    return std::ceil(std::log10(value)/logG);
}

double getCurrentAlpha(double alpha) {
    return 2*(alpha/(1 + pow(alpha,2)));
}

double getCurrentGamma(double alpha) {
    return (1+alpha)/(1-alpha);
}

double getCurrentLogG(double gamma) {
    return std::log10(gamma);
}



double getQuantileFraction(int kth, int I) {

    double q = (double)floor( (kth-1)*100/(I-1) );     
    q /= 100.0;
    if (q < 0.0 || q > 1.0) {
        std::cerr << "ERROR while computing q value for DDSketch\n";
        exit(1);
    }
    return q;
}


double estimateQ(std::map<int, int>& Sketch, double q, double gamma, int n) {

    double estimate = 0.0;
    double fraction = q*(n-1);
    
    
    std::map<int, int>::iterator it = Sketch.begin();
    int i = it->first;
    int count = it->second;
    
    while (count <= fraction) {
        ++it;
        i = it->first;
        count += it->second;
    }//wend
    
    estimate = (2.0 * pow(gamma,i))/(gamma+1.0);
    
    #ifdef DEBUG
        std::cout << "Quantile: " << estimate << ", count = "<<count << ", i = "<< i << ", population " << n << ", fraction " << fraction<< std::endl;
    #endif

    return estimate; 
}


void debugSketch(std::map< int, int>& mySketch) {

    fprintf(stdout,"\nSketch is : \n\t Key \t Count\n");
    int TotalCount = accumulate( mySketch.begin(), mySketch.end(), 0, []( int acc, std::pair<int, int> p ) { return ( acc + p.second ); } );
    
    int loop= 1;
    for(auto it=mySketch.begin(); it != mySketch.end(); ++it) {
        fprintf(stdout,"%d) \t%+12d, \t%d\n", loop++, it->first, it->second);
    }//for
    fprintf(stdout,"Total differences contained in sketch %d, over %lu buckets\n\n", TotalCount, mySketch.size());
}




double estimator(std::map< int, int>&mySketch, double q, double gamma) {

    double sum = 0;
    std::map< int, int>::iterator it = mySketch.begin();

    int i = it->first;
    double b_count = it->second;
    sum = b_count;
    while(sum<=q && it != mySketch.end()) {
        ++it;

        i = it->first;
        sum += it->second;
    }//wend

    //Return the estimation x_q of bucket with key i
    return (2.0 * pow(gamma,i))/(gamma+1.0);
}

void logQuantiles(FILE *fp, std::map< int,int>&mySketch, int collapses, double gamma, double *exactDiffs, int len) {
    
    double I = 1.0*len; 
    double population = 0;
    std::map< int,int >::iterator it;
    for (it = mySketch.begin(); it != mySketch.end(); ++it) {
        population +=it->second;
    }//for (it)

    double Amin = estimator(mySketch, 0.00,  gamma);
    double Aq1 = estimator(mySketch, 0.25*(population-1),  gamma);
    double Aq2 = estimator(mySketch, 0.50*(population-1),  gamma);
    double Aq3 = estimator(mySketch, 0.75*(population-1),  gamma);
    double Amax = estimator(mySketch, 1.00*(population-1),  gamma);

    double Emin = quickselect(exactDiffs, len, 0); 
    double Eq1 = quickselect(exactDiffs, len, (int)std::ceil(I/4.0)-1 ); 
    double Eq2 = quickselect(exactDiffs, len, (int)std::ceil(I/2.0)-1 ); 
    double Eq3 = quickselect(exactDiffs, len, (int)std::ceil(3*I/4.0)-1 ); 
    double Emax = quickselect(exactDiffs, len, len-1); 

    fprintf(fp, "%.0f,%lu,%d,", population, mySketch.size(), collapses);
    fprintf(fp, "%.6f,%.6f,%.6f,%d,", Emin, Amin, std::abs((Amin-Emin)/Emin),0);    
    fprintf(fp, "%.6f,%.6f,%.6f,%d,", Eq1, Aq1, std::abs((Aq1-Eq1)/Eq1),(int)std::ceil(I/4.0)-1);
    fprintf(fp, "%.6f,%.6f,%.6f,%d,", Eq2, Aq2, std::abs((Aq2-Eq2)/Eq2),(int)std::ceil(I/2.0)-1);
    fprintf(fp, "%.6f,%.6f,%.6f,%d,", Eq3, Aq3, std::abs((Aq3-Eq3)/Eq3),(int)std::ceil(3*I/4.0)-1);
    fprintf(fp, "%.6f,%.6f,%.6f,%d\n", Emax, Amax, std::abs((Amax-Emax)/Emax),len-1);
}



//****** ****** ****** ****** ****** ************ ************ ************ ************ ****** Uniform Collapse of the sketch

void uniformCollapse(std::map<int, int>& mySketch) {
 
    std::map<int, int> newSketch; 
    
    std::map<int, int>::iterator it = mySketch.begin();
    while (it != mySketch.end() ) {

        int key = it->first;
        std::map<int, int>::iterator it_collapsed;

        if (key%2 != 0) 
        {
            
            int newK = (key+1)/2;

            it_collapsed = mySketch.find(key+1);
            if (it_collapsed == mySketch.end()) { 
                
                newSketch[newK] = it->second;
               
                it = mySketch.erase(it);
            } else {
                newSketch[newK] = it->second + it_collapsed->second; 
                it = mySketch.erase(it); 
                it = mySketch.erase(it); 

            }// fi it_collapsed

        }//fi key odd 
        else 
        {   
            int newK = (key)/2; 
            
            newSketch[newK] = it->second; 
            
            it = mySketch.erase(it); 
        }//fi key even

    }//wend

    mySketch.swap(newSketch);
}


void collapseUniformly(std::map<int, int>& mySketch) {

    std::map<int, int> newSketch; 
    std::map<int, int>::iterator it;

    it = mySketch.begin();
    if (it->first == -MIN_KEY) { 
        newSketch[-MIN_KEY] += it->second;
        ++it;
    }

    while(it != mySketch.end()){
        
        double k = it->first;
        int k_new = std::ceil(k/2);
        newSketch[k_new] += it->second;
        
        ++it;
    }//wend

    mySketch.swap(newSketch);
}


int performCollapse(std::map<int, int>& Sketch, int sketchBound, double *currentAlpha, double *currentGamma, double *currentLogG, int *SketchSize) {

    int collapse_executed = 0;
    int currentSize = Sketch.size();

    while (currentSize > sketchBound) { 
    
        *currentAlpha = getCurrentAlpha(*currentAlpha);
        *currentGamma = getCurrentGamma(*currentAlpha);
        *currentLogG  = getCurrentLogG(*currentGamma);

        collapseUniformly(Sketch); 
        ++collapse_executed; 
        
        currentSize = Sketch.size();

        #ifdef DEBUG
            std::cout << "New error params [current (α, 𝛄, LogG) " << *currentAlpha;
            std::cout << ", " << *currentGamma;
            std::cout << ", " << *currentLogG << " ]" << std::endl;
        
            int TotalCount = accumulate( Sketch.begin(), Sketch.end(), 0, []( int acc, std::pair<int, int> p ) { return ( acc + p.second ); } );
            std::cout<< "After " << collapse_executed << " collapse() the sketch size is "<< Sketch.size();
            std::cout<< "\tTotal count: " << TotalCount << std::endl;
        #endif

    }//wend - collapse()
    
    (*SketchSize) = currentSize;
    return collapse_executed;
}


//****** ****** ****** ****** ****** ************ ************ ************ ************ ****** Filling the sketch

int fillSketch(int pos, double *window, double gamma, double LogG, std::map<int, int>& Sketch) {

    double new_diff = -1;
    int currentBi;
    int count_added = 0;

   
    for(int j = pos-1; j>=0; --j) {            

        new_diff = std::abs(window[pos]-window[j]);

        currentBi = getKeyFor(new_diff, gamma, LogG);   

        Sketch[currentBi] += 1;

        ++count_added;                                              
    }//for ADD diffs

    return count_added;
}


//****** ****** ****** ****** ****** ************ ************ ************ ************ ****** Update the sketch

int uniformRemove(std::map<int,int> &sketch, double gamma, double logGamma, int pos, int ndiffs, double *Pwindow, int s) {
    
    double old_item = Pwindow[pos];
    int sample = (s-1)/ndiffs; //sampling frequency: one difference every "sample" items
    int key,res;

    int r = pos + 1;
    int l = pos - 1;
    int count = 0, i=0;
    while ((count < ndiffs) && (l>=0 || r<s)) {    

            while ((r < s) && (count<ndiffs)) {

                key = getKeyFor(std::abs(Pwindow[r]-old_item), gamma, logGamma);
                res = decreaseBinCount(key, sketch);

                if (res == 1) {
                    ++count;
                }//fi res 

                r += sample;
            }//wend r
           
            while ((l >= 0) && (count<ndiffs)) {

                key = getKeyFor(std::abs(Pwindow[l]-old_item), gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                
                if (res == 1) {
                    ++count;
                }//fi res 
                
                l -= sample;
            }//wend l

        ++i;
        r = pos + 1 + i;
        l = pos - 1 - i;
    }//wend
    
    return count;
}

int uniformAdd(std::map<int,int> &sketch, double gamma, double logGamma, int pos, int ndiffs, double *Pwindow, int s) {

    double new_item = Pwindow[pos];
    
    int r = pos + 1;
    int l = pos - 1;
    int count = 0;
    int key;

    int sample = (s-1)/ndiffs; 
            
            while (r < s && count<ndiffs){
                key = getKeyFor(std::abs(Pwindow[r]-new_item), gamma, logGamma);
                sketch[key] += 1;
                
                r+=sample;
                ++count;
                
            }
           
            while (l >= 0 && count<ndiffs){
                key = getKeyFor(std::abs(Pwindow[l]-new_item), gamma, logGamma);
                sketch[key] += 1;
                
                l-=sample;
                ++count;
            
            }//wend l

    return count;
}



int decreaseBinCount(int key, std::map<int, int>& sketch) {
    
    int res = 1;
    std::map<int, int>::iterator it = sketch.find(key);

    if ( it == sketch.end() ) {
        res = -1;
        #ifdef PARTIAL
            return res;
        #else // ALL (s-1) differences
            std::cout<<"removeWeightToSketch(): ERROR: try to remove item from non existing bucket "<< key << "\n";
            exit(1);
        #endif
    }
    
    it->second -= 1;
    if (it->second == 0) {
        sketch.erase(it);
    }
    return res;
}



int selectDiffsToRemove2(std::map<int,int> &sketch, double gamma, double logGamma, int pos, int ndiffs, double *Pwindow, int s) {
    
    double old_item = Pwindow[pos];
    int key, res;
    
    int r = pos + 1;
    int l = pos - 1;
    int count = 0, i=0;
    int quit = 0;
    
   
    while (count < ndiffs && !quit) {

        if (r<s && l>=0) {

            
            double d1 = std::abs(Pwindow[l] - old_item);
            double d2 = std::abs(old_item - Pwindow[r]);

            if ( d1 <= d2 ){
                key = getKeyFor(d1, gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                if (res == 1) {
                    ++count;    
                }
                --l;
            } else {
                key = getKeyFor(d2, gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                if (res == 1) {
                    ++count;
                }
                ++r;
            }//fi smallest diff
        } else {
            quit = 1;
            while ((r<s) && (count<ndiffs)) {
                
                key = getKeyFor(std::abs(Pwindow[r]-old_item), gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                if (res == 1) {
                    ++count;
                }//fi res
                ++r;
            }//wend r

            while ((l>=0) && (count<ndiffs)) {

                key = getKeyFor(std::abs(Pwindow[l]-old_item), gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                if (res == 1) {
                    ++count;
                }//fi res
                --l;
            }//wend l

        }//fi (r && l)
    }//wend
    
    return count;
}



int selectDiffsToRemove(std::map<int,int> &sketch, double gamma, double logGamma, int pos, int ndiffs, double *Pwindow, int s) {
    
    double old_item = Pwindow[pos];
    int key, res;
    
    int r = pos + 1;
    int l = pos - 1;
    int count = 0, i=0;

    while ((count < ndiffs) && (l>=0 || r<s)) {
        
        #ifdef DEBUG
            std::cout <<"Difference(-) #"<< count+1 << "): " << Pwindow[l] << " | " << old_item << " | " << Pwindow[r];
        #endif

        if ((r<s && l>=0) && ( (std::abs(Pwindow[l]-old_item )) <= std::abs(old_item - Pwindow[r]) ) ) {

            #ifdef DEBUG
            std::cout << ": choose left " << std::endl;
            #endif

            key = getKeyFor(std::abs(old_item-Pwindow[l]), gamma, logGamma);
            res = decreaseBinCount(key, sketch);
            if (res == 1) {
                ++count;    //PARTIAL mode
            }
            --l;
            
        } else if ((r<s && l>=0) && ( (std::abs(Pwindow[l]-old_item )) > std::abs(old_item - Pwindow[r]) ) ) {

            #ifdef DEBUG
            std::cout << ": choose right " << std::endl;
            #endif
            
            key = getKeyFor(std::abs(Pwindow[r]-old_item), gamma, logGamma);
            res = decreaseBinCount(key, sketch);
            if (res == 1) {
                ++count;
            }
            ++r;
        
        } 
        else if (l<0) {
            
            #ifdef DEBUG
                std::cout << ": choose right from now on" << std::endl;
            #endif

            while ((r<s) && (count<ndiffs)) {

                key = getKeyFor(std::abs(Pwindow[r]-old_item), gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                if (res == 1) {
                    ++count;
                }//fi res
                ++r;
            }//wend r
        
        } else if (r>=s) {
           
            #ifdef DEBUG
                std::cout << ": choose left from now on" << std::endl;
            #endif


            while ((l>=0) && (count<ndiffs)) {

                key = getKeyFor(std::abs(Pwindow[l]-old_item), gamma, logGamma);
                res = decreaseBinCount(key, sketch);
                if (res == 1) {
                    ++count;
                }//fi res
                --l;
            }//wend l

        }//fi cases

    }//wend
    
    return count;
}



int selectDiffsToAdd(std::map<int,int> &sketch, double gamma, double logGamma, int pos, int ndiffs, double *Pwindow, int s) {

    double new_item = Pwindow[pos];
    
    int r = pos + 1;
    int l = pos - 1;
    
    int count = 0;
    int key;
    while ((count < ndiffs) && (l>=0 || r<s)) 
    { 
        #ifdef DEBUG
        std::cout <<"Difference(+) "<< count+1 << "): " << Pwindow[l] << " | " << new_item << " | " << Pwindow[r];
        #endif

        if ((r<s && l>=0) && ( (std::abs(Pwindow[l]-new_item )) <= std::abs(new_item - Pwindow[r]) ) ) {
            
            #ifdef DEBUG
            std::cout << ": choose left" << std::endl;
            #endif
            
            key = getKeyFor(std::abs(new_item-Pwindow[l]), gamma, logGamma);
            sketch[key] += 1;
            ++count;
            --l;
        }
        else if ((r<s && l>=0) && ( (std::abs(Pwindow[l]-new_item )) > std::abs(new_item - Pwindow[r]) ) ) {
            
            #ifdef DEBUG
            std::cout << ": choose right" << std::endl;
            #endif

            key = getKeyFor(std::abs(Pwindow[r]-new_item), gamma, logGamma);
            sketch[key] += 1;
            ++count;
            ++r;
        } 
        
        else if (l<0) {
            #ifdef DEBUG
            std::cout << ": choose right from now on" << std::endl;
            #endif
    

            for(int i = r; i < s; ++i) {

                key = getKeyFor(std::abs(Pwindow[i]-new_item), gamma, logGamma);
                sketch[key] += 1;
                ++count;
                if (count == ndiffs){
                    return count;
                }
            }//for 
            r = s;
        } 
        else if (r>=s) {
           
            #ifdef DEBUG
            std::cout << ": choose left from now on" << std::endl;
            #endif

           for(int i = l; i >= 0; --i) {

                key = getKeyFor(std::abs(new_item-Pwindow[i]), gamma, logGamma);
                sketch[key] += 1;
                ++count;
                if (count == ndiffs){
                    return count;
                }
            }//for
            l = -1;

        }//fi cases

    }//wend
    return count;
}



int selectDiffsToAdd2(std::map<int,int> &sketch, double gamma, double logGamma, int pos, int ndiffs, double *Pwindow, int s) {

    double new_item = Pwindow[pos];
    
    int r = pos + 1;
    int l = pos - 1;
    int count = 0;
    int key;

    while (count < ndiffs) { 
        
        if (r<s && l>=0) {
            double d1 = std::abs(Pwindow[l]-new_item); 
            double d2 = std::abs(new_item - Pwindow[r]);
            
            if (d1<=d2){
                key = getKeyFor(d1, gamma, logGamma);
                sketch[key] += 1;
                ++count;
                --l;
            } else {
                key = getKeyFor(d2, gamma, logGamma);
                sketch[key] += 1;
                ++count;
                ++r;
            }//fi smallest diff
        } else {
            while (r<s && count < ndiffs){
                key = getKeyFor(std::abs(Pwindow[r]-new_item), gamma, logGamma);
                sketch[key] += 1;
                ++count;
                ++r;
            }//wend r

            while (l>=0 && count < ndiffs){
                key = getKeyFor(std::abs(new_item-Pwindow[l]), gamma, logGamma);
                sketch[key] += 1;
                ++count;
                --l;
            }//wend l

        }//fi (r && l)

    }//wend
    return count;
}





int updateSketch(double old_item, double new_item, double *Pwindow, int s, std::map<int,int>& Sketch, double gamma, double logGamma, int ndiffs) {
    
    int pos = -1;               
    int population = 0;         
    int removed = 0, added = 0;


    // 1. equal items
    if (old_item == new_item) {
        return population;
    }//fi old=new


    if (old_item < new_item) {

        pos = bsearch(old_item, Pwindow, s);
        if (pos == -1){
            std::cerr << "ERROR while searching an existing item in Pwindow " << std::endl;
            exit(1);
        }

        #ifndef UNIFORM 
            removed = selectDiffsToRemove2(Sketch, gamma, logGamma, pos, ndiffs, Pwindow, s);
        #else
            removed = uniformRemove(Sketch, gamma, logGamma, pos, ndiffs, Pwindow, s);
        #endif
       
        for ( int p = pos; p < s-1; ++p) {

            if (Pwindow[p+1] < new_item) {
                Pwindow[p] = Pwindow[p+1];
            } else {
                Pwindow[p] = new_item;
                
                #ifndef UNIFORM 
                    added = selectDiffsToAdd2(Sketch, gamma, logGamma, p, removed, Pwindow, s);
                #else
                    added = uniformAdd(Sketch, gamma, logGamma, p, removed, Pwindow, s);
                #endif

                return (added-removed); 
            }//fi check

        }//for insert new_item

        Pwindow[s-1] = new_item;

        #ifndef UNIFORM 
            added = selectDiffsToAdd2(Sketch, gamma, logGamma, s-1, removed, Pwindow, s);
        #else
            added = uniformAdd(Sketch, gamma, logGamma, s-1, removed, Pwindow, s);
        #endif

        return (added-removed); 
    }// fi (old_item < new_item)
    else 
    { 
        pos = bsearch(old_item, Pwindow, s);
        if (pos == -1){
            std::cerr << "ERROR while searching an existing item in Pwindow " << std::endl;
            exit(1);
        }

        
        #ifndef UNIFORM 
           removed = selectDiffsToRemove2(Sketch, gamma, logGamma, pos, ndiffs, Pwindow, s);
        #else
            removed = uniformRemove(Sketch, gamma, logGamma, pos, ndiffs, Pwindow, s);
        #endif

        for (int p = pos; p > 0; --p) {

            if (Pwindow[p-1] > new_item) {
                Pwindow[p] = Pwindow[p-1];
            } else {
                
                Pwindow[p] = new_item;
                
                #ifndef UNIFORM 
                    added = selectDiffsToAdd2(Sketch, gamma, logGamma, p, removed, Pwindow, s);
                #else
                    added = uniformAdd(Sketch, gamma, logGamma, p, removed, Pwindow, s);
                #endif

                return (added-removed); 
            }//fi check
        }//for

        Pwindow[0] = new_item;

        #ifndef UNIFORM 
            added = selectDiffsToAdd2(Sketch, gamma, logGamma, 0, removed, Pwindow, s);
        #else
            added = uniformAdd(Sketch, gamma, logGamma, 0, removed, Pwindow, s);
        #endif

        return (added-removed); 
    }// fi (new_item < old_item)  

    return population;
} 

//********************************************************************************************

int computeDiffs(double new_item, double old_item, double Pitem, double gamma, double logG, std::map<int,int>&sketch) {
    
    int added = 0;
    int removed = 0;

    double diffA = std::abs(Pitem - new_item); 
    int keyA = getKeyFor(diffA, gamma, logG);
                
    double diffR = std::abs(Pitem - old_item);
    int keyR = getKeyFor(diffR, gamma, logG);
            
    if (keyA != keyR) {
        sketch[keyA] += 1;
        ++added;

        std::map<int,int>::iterator it = sketch.find(keyR);
        if (it != sketch.end()) {
            ++removed;
            it->second -= 1;
            if (!it->second){
                sketch.erase(it);
            }
        } else {
            std::cerr << "ERROR : key not found in sketch while deleting oldest item"<<std::endl;
            exit(1);
        }
    }//keyA!=keyR

    return added-removed;
}



void updateSynopsis(double old_item, double new_item, double *Pwindow, int s, std::map<int,int>& Sketch, double gamma, double logGamma) {
    
    int posA = -1;                       
    int posR = -1;                       
    int removed = 0, added = 0;          
    bool deleted = false;                
    bool inserted = false;               
    
    int balance = 0;
    int counter = 0;

    if (old_item < new_item) {
        int p = 0;              
        
        while (!deleted) {
            
            if (Pwindow[p] == old_item) 
            {
                posR = p; 

                if ( ((p<s-1) && (new_item <= Pwindow[p+1]) ) || (p == s-1) ) {
                    posA = p;
                    Pwindow[p] = new_item;
                    inserted = true;
                } else {

                    Pwindow[p] = Pwindow[p+1];
                    
                    balance = computeDiffs(new_item, old_item, Pwindow[p], gamma, logGamma, Sketch);
                    if (balance) {
                        std::cerr<<"ERROR on updating the sketch\n";
                        exit(1);
                    }
                    ++counter;
                }//fi new_item
                deleted = true;
            } 
            else 
            {
                balance = computeDiffs(new_item, old_item, Pwindow[p], gamma, logGamma, Sketch);
                if (balance) {
                    std::cerr<<"ERROR on updating the sketch\n";
                    exit(1);
                }
                ++counter;
            }//fi
            ++p;
        }//wend

        while (!inserted) {

            if ( ((p<s-1) && (new_item <= Pwindow[p+1]) ) || (p == s-1) ) {
                posA = p;
                Pwindow[p] = new_item;
                inserted = true;
            } else {

                Pwindow[p] = Pwindow[p+1];
                
                balance = computeDiffs(new_item, old_item, Pwindow[p], gamma, logGamma, Sketch);
                if (balance){
                    std::cerr<<"ERROR on updating the sketch\n";
                    exit(1);
                }
                ++counter;
            }//fi new_item
            ++p;
        }

        while (p<s) {
            balance = computeDiffs(new_item, old_item, Pwindow[p], gamma, logGamma, Sketch);
            if (balance){
                std::cerr<<"ERROR on updating the sketch\n";
                exit(1);
            }
            ++counter;
            ++p;
        }//wend

    } else {
        int q = s-1; 

        while (!deleted) {
            
            if (Pwindow[q] == old_item) 
            {
                posR = q;

                if ( ((q>0) && (new_item >= Pwindow[q-1]) ) || (q == 0) ) {
                    posA = q;
                    Pwindow[q] = new_item;
                    inserted = true;
                } else {

                    Pwindow[q] = Pwindow[q-1];
                    
                    balance = computeDiffs(new_item, old_item, Pwindow[q], gamma, logGamma, Sketch);
                    if (balance){
                        std::cerr<<"ERROR on updating the sketch\n";
                        exit(1);
                    }
                    ++counter;
                }//fi new_item
                deleted = true;
            } 
            else 
            {
                balance = computeDiffs(new_item, old_item, Pwindow[q], gamma, logGamma, Sketch);
                if (balance) {
                    std::cerr<<"ERROR on updating the sketch\n";
                    exit(1);
                }
                ++counter;
            }//fi
            --q;
        }//wend
        
        while (!inserted) {
            if ( ((q>0) && (new_item >= Pwindow[q-1]) ) || (q == 0) ) {
                posA = q;
                Pwindow[q] = new_item;
                inserted = true;
            } else {

                Pwindow[q] = Pwindow[q-1];
                
                balance = computeDiffs(new_item, old_item, Pwindow[q], gamma, logGamma, Sketch);
                if (balance){
                    std::cerr<<"ERROR on updating the sketch\n";
                    exit(1);
                }
                ++counter;
            }//fi new_item
            --q;
        }

        while (q>-1) {
            balance = computeDiffs(new_item, old_item, Pwindow[q], gamma, logGamma, Sketch);
            if (balance){
                std::cerr<<"ERROR on updating the sketch\n";
                exit(1);
            }
            ++counter;
            --q;
        }//wend


    }//fi 
} 
