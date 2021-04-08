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

#include "IIS.h"
#include <iostream>

//********************************************************************************** INSERTION SORT V5



int isort_v5(double *V, int len, double new_item) {
    
    if (len == 0) {
        V[0] = new_item;
        return 0; 
    }

    double tmp = 0;
    int k = 0;
    
    int pos = -1;               
    
    for(k = 0; k < len; ++k){

        if (new_item <= V[k]){

            if (pos == -1) {
                pos = k;
            } 
            
            tmp = V[k];  
            V[k] = new_item;
            new_item = tmp;          
        }

    }//for
    
    if (pos == -1) {
        pos = k;
    }
    V[k] = new_item;

    return pos;
}






int updateValues_v5(double *V, int n, double new_item, double old_item, Pos *positions){
    int pos = -1;              
    if (old_item == new_item)
    {
         for(int k = 0; k < n; ++k){
            if (old_item == V[k]){
                pos = k;
                break;          
            }
        }//for
        positions->_newpos = pos;
        positions->_oldpos = pos;
        return pos;
    } 

    if (new_item < old_item) {

        double tmp = old_item;        
        for(int k = 0; k < n; ++k) {

            if (new_item <= V[k]) {
                
                if (pos == -1){
                    pos = k;
                    positions->_newpos = pos;
                }
                
                tmp = V[k];  
                V[k] = new_item;
                new_item = tmp;

                if (new_item==old_item){
                    positions->_oldpos = k;
                    break;
                }
            }//fi

        }//for k
    
    } else {
        double tmp;
        for(int k = n-1; k >= 0; --k){
            if (new_item >= V[k]) {
                
                if (pos == -1){
                    pos = k;
                    positions->_newpos = pos;
                }

                tmp = V[k];
                V[k] = new_item;
                new_item = tmp;
            
                if (new_item == old_item){
                    positions->_oldpos = k;
                    break;
                }
            }//fi
        }//for k
    }//fi old_item

return pos;
}



int bsearch(double value, double *p, int size) {
    
    int l = 0;
    int r = size-1;
    while (l <= r) { 
        int m = l + (r - l) / 2; 
  
        if (p[m] == value) 
            return m; 
  
        if (p[m] < value) {
            l = m + 1; 
        }
        else {
            r = m - 1; 
        }
    }//wend 
    return -1;
}



void updateSortedWindow(double *Pi, int s, double new_item, double old_item) {
    
    int pos = -1;

    if (old_item == new_item) {
        return;
    }//fi old=new


    if (old_item < new_item) {   
        pos = bsearch(old_item, Pi, s);
        if (pos == -1){
            std::cerr << "ERROR while searching an existing item in Pwindow " << std::endl;
            exit(1);
        }
       
        for ( int p = pos; p < s-1; ++p) {
            if (Pi[p+1] < new_item){
                Pi[p] = Pi[p+1];
            } else {
                Pi[p] = new_item;
                return; 
            }//fi check
        }//for insert new_item

        Pi[s-1] = new_item;
        return;
    }// fi (old_item < new_item)
    else 
    { 
        pos = bsearch(old_item, Pi, s);
        if (pos == -1){
            std::cerr << "ERROR while searching an existing item in Pi " << std::endl;
            exit(1);
        }

        for (int p = pos; p > 0; --p) {
            if (Pi[p-1] > new_item){
                Pi[p] = Pi[p-1];
            } else {
                Pi[p] = new_item;
                return; 
            }//fi check
        }//for

        Pi[0] = new_item;

        return;
    }// fi (new_item < old_item)  
} 
