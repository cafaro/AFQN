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


#ifndef __IIS_H__
#define __IIS_H__

typedef struct Pos {

    int _oldpos;
    int _newpos;
}Pos; 




int isort_v5(double *V, int len, double new_item);


int updateValues_v5(double *V, int n, double new_item, double old_item, Pos *positions);



// ******************* Working on sorted permutation of W
int bsearch(double value, double *p, int size);


void updateSortedWindow(double *Pi, int s, double new_item, double old_item);


#endif //__IIS_H__
