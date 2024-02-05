#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

//copy form sa_moral_0326R3.cpp on April 8 2016 in Madison
extern "C" SEXP  Structural_Hamming_Distance(SEXP H_, SEXP G_, SEXP p_){
/*%Input:
%H is the PDAG by learning algorithm
%G is the true PDAG representing equivalence

%Output:
%miss: the number of missing edges in H
%extra: the number of extra edges in H
%SHD: the Structural Hamming Distance  defined in page 30, i.e., Algorithm 4: SHD Algorithm, by 
%Tsamardinos, I., Brown, L.E., Aliferis, C.F.: 
%The max-min hill-climbing Bayesian network structure learning algorithm. Machine Learning 65(1), 31o?=C78 (2006)*/
  double * G = REAL(G_);
  double * H = REAL(H_);
  int p = INTEGER(p_)[0];// the number of vairables
  double *diff_graph = new double[p*p];
  //diff_graph = abs( G - H );
  for(int i=0; i<p*p; i++){
    diff_graph[i] = fabs(G[i] - H[i]);
  }
  
  
  double me = 0; //miss edges
  double ee = 0; //extra edges
  double md = 0; //miss direction
  double ed = 0; //extra direction
  double rd = 0; //reverse direction
  for(int i=0; i<p; i++){
    for(int j=0; j<p; j++){
      if(diff_graph[i+j*p]!=0){
        if(H[i+j*p] == 0 && H[j+i*p] == 0){
            if(diff_graph[i+j*p] != diff_graph[j+i*p])
                me = me + 1;
            else
                me = me + 0.5;            
        }            
        else if(G[i+j*p] == 0 && G[j+i*p] == 0){
            if(diff_graph[i+j*p] != diff_graph[j+i*p])
                ee = ee + 1;
            else
                ee = ee + 0.5;            
        }            
        else{
            if(G[i+j*p]==1 && G[j+i*p] ==1){
                ed = ed + 1;
            }                
            else if(G[i+j*p] ==1 && G[j+i*p] ==0){
                if( H[i+j*p] ==1 && H[j+i*p]==1){
                    md = md + 1;
                }                    
                else if( H[i+j*p] ==0 && H[j+i*p]==1){
                    rd = rd + 1;
                }
            }                
            else if(G[i+j*p] ==0 && G[j+i*p] ==1){
                if(H[i+j*p] ==1 && H[j+i*p]==1){
                    md = md + 1;
                }                    
                else if( H[i+j*p] ==1 && H[j+i*p]==0 ){
                    rd = rd + 1;
                }
            }
        }
      }        
    }
  }
  
  delete []diff_graph;
  
  double SHD = ee + me  + md + ed + rd;
  //diff = [SHD, me, ee, md, ed, rd];
  
  
  SEXP diff;
  PROTECT(diff=allocVector(REALSXP, 6));
  double *rdiff=REAL(diff);
  rdiff[0] = SHD;
  rdiff[1] = me;
  rdiff[2] = ee;
  rdiff[3] = md;
  rdiff[4] = ed;
  rdiff[5] = rd;
  
  UNPROTECT(1);
    return diff;
}


// 
// shd_cpdag<-function(G, H){
// #H is the PDAG by learning algorithm
// #G is the true PDAG representing equivalence
//   p = ncol(G);
//   if(ncol(G)!=ncol(H) ||nrow(G)!=nrow(H) ||nrow(G)!=ncol(H) )
//     stop("G and H must be two square matrix with the same number of rows\n");
//   dd<-.Call("Structural_Hamming_Distance", as.double(H), as.double(G), as.integer(p) );
//   d<-matrix(dd, ncol=6, byrow=TRUE);
//   return(d)
// }