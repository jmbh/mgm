#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix mMRFCsampler(NumericMatrix Data, int n, int nNodes, NumericVector type_c,
                           NumericVector levels, int nIter, NumericMatrix thresh_m,
                           NumericMatrix graphe, IntegerVector inde) {
  
  NumericVector indicatorstate;
  NumericVector potcat;
  NumericVector potcon;
  
  for(int p=0; p<n; p++) { // loop cases
    
    
    // generate initial states
    for (int node=0; node<nNodes; node++) {
      if(type_c[node] == 1)
        Data(p,node) =  (double)round(R::runif(0.5, levels[node]+.5));
      if(type_c[node] == 2)
        Data(p,node) =  (double)R::rnorm(1,1);
      if(type_c[node] == 3)
        Data(p,node) =  (double)R::rpois(1);
      if(type_c[node] == 4)
        Data(p,node) =  (double)R::rexp(1);
    }
    
    
    for(int iter = 0; iter<nIter; iter++) { // loop iterations
      
      for(int node=0; node<nNodes; node++) { // loop nodes
        
        // CATEGORICAL
        if(type_c[node]==1)
        {
          
          NumericVector potential_stor(levels[node]); //storage for potential of each category
          
          for(int l = 0; l<levels[node]; l++) { // loop categories
            
            
            //create storage for interaction terms
            // gives an indicator how many cat vs. cont variables (to create approporate storage)
            NumericVector ind_type(nNodes);
            for(int i=0; i<nNodes; i++)
            {
              if(type_c[i]==1)
                ind_type[i] = 1;
              else
                ind_type[i] = 0;
            }
            //create storage objects
            potcat = rep(0,(sum(ind_type)-1)); //we know we deal type_c[node] = categorical
            potcon = rep(0,(nNodes - sum(ind_type)));
            int cat = 0; //"manual" indices, because depends on sequence of con-cat-....
            int con = 0;
            
            for(int node2 = 0; node2<nNodes; node2++) { // loop other nodes
              
              if(node2!=node) {
                
                if(type_c[node2]==1)
                {
                  
                  // which categoryy is present
                  indicatorstate = rep(0,levels[node2]);
                  for(int k=0;k<levels[node2]; k++) {
                    if(Data(p,node2)==(k+1)) indicatorstate[k]=1;
                  }
                  
                  
                  // cut out correct part of graphe
                  IntegerVector ind_rows(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}
                  
                  IntegerVector ind_cols(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}
                  
                  NumericMatrix n_m1(sum(ind_rows), ind_cols.size());
                  
                  int rowc = 0;
                  for(int r = 0; r<ind_cols.size(); r++) {
                    if(ind_rows[r]==1) {
                      n_m1(rowc,_) = graphe(r,_);
                      rowc++;
                    }
                  }
                  
                  NumericVector n_v1 = n_m1(l,_); //cut out level/state
                  NumericVector n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)
                  
                  potcat[cat] = sum(n_v2 * indicatorstate); // interaction term of node2 categorical variable
                  cat++;
                  
                  
                } else { //if node2=(con)
                  
                  // cut out correct part of graphe (same as above)
                  IntegerVector ind_rows(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}
                  
                  IntegerVector ind_cols(inde.size());
                  for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}
                  
                  NumericMatrix n_m1(sum(ind_rows), ind_cols.size());
                  
                  int rowc = 0;
                  for(int r = 0; r<ind_cols.size(); r++) {
                    if(ind_rows[r]==1) {
                      n_m1(rowc,_) = graphe(r,_);
                      rowc++;
                    }
                  }
                  
                  NumericVector n_v1 = n_m1(l,_); //cut out level/state
                  
                  NumericVector n_v2;
                  n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)
                  
                  
                  NumericVector dummy_postcon; //apparently necessary to get from vector->double
                  dummy_postcon = n_v2 * Data(p,node2); // interaction term of node2 (continuous variable)
                  potcon[con] = dummy_postcon[0]; // interaction term of node2 (continuous variable)
                  con++;
                  
                } // end: if node2: cond vs. con
                
              } // end: if node2!=node
              
            } // end: other nodes
            
            double potraw =  thresh_m(node,l) + sum(potcat) + sum(potcon);
            
            // to avoid infinite values when taking exp(potraw); this unelegant thresholdung makes no difference in our setting
            if(potraw>100) {potraw = 100;}
            
            //debug
           // Rcpp::Rcout << "iter " << iter << " case " << p << " node " << node << " potraw " << potraw << std::endl;
            
            potential_stor[l] =  exp(potraw);
            
          } // end: loop categories
          
          //sample state proportional to potentials
          double samplesum = sum(potential_stor);
          
          
          double randomNumber = R::runif(0, samplesum);
          double newsum = 0;
          for (int ind=0; ; ind++) {
            newsum += potential_stor[ind];
            if (randomNumber <= newsum) {
              Data(p,node) = ind+1;
              break;
            }
          }
          
        } else { // CONTINUOUS (most code just copied from above)
          
          
          //create storage for interaction terms (this could be outside)
          // gives an indicator how many cat vs. cont variables (to create approporate storage)
          NumericVector ind_type(nNodes);
          for(int i=0; i<nNodes; i++)
          {
            if(type_c[i]==1)
              ind_type[i] = 1;
            else
              ind_type[i] = 0;
          }
          //create storage objects
          potcat = rep(0,(sum(ind_type)));
          potcon = rep(0,(nNodes - sum(ind_type)-1)); //we know we deal type_c[node] = continuous
          int cat = 0; //"manual" indices, because depends on sequence of con-cat-....
          int con = 0;
          
          for(int node2 = 0; node2<nNodes; node2++) { // loop other nodes
            
            if(node2!=node) {
              
              //node2 = cat
              if(type_c[node2]==1) {
                
                
                // which categoryy is present
                indicatorstate = rep(0,levels[node2]);
                for(int k=0;k<levels[node2]; k++) {
                  if(Data(p,node2)==(k+1)) indicatorstate[k]=1;
                }
                
                
                // cut out correct part of graphe
                IntegerVector ind_rows(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;}
                
                IntegerVector ind_cols(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;}
                
                NumericMatrix n_m1(sum(ind_rows), ind_cols.size());
                
                int rowc = 0;
                for(int r = 0; r<ind_cols.size(); r++) {
                  if(ind_rows[r]==1) {
                    n_m1(rowc,_) = graphe(r,_);
                    rowc++;
                  }
                }
                
                NumericVector n_v1 = n_m1(0,_); //to use the same structure is above; there is of course only 1 row
                NumericVector n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2)
                
                potcat[cat] = sum(n_v2 * indicatorstate); // interaction term of node2 categorical variable
                cat++;
                
                
              } else { //node2 = con
                
                // cut out correct part of graphe (same as above)
                IntegerVector ind_rows(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node+1)) ind_rows[k]=1;} // present node is flagged with 1 (1 integer in this case)
                
                IntegerVector ind_cols(inde.size());
                for(int k=0;k<inde.size(); k++) {if(inde[k]==(node2+1)) ind_cols[k]=1;} // second node2 is flagged with 1 (1 integer in this case)
                
                NumericMatrix n_m1(sum(ind_rows), ind_cols.size()); // in this case always a 1x1 matrix
                
                
                int rowc = 0;
                for(int r = 0; r<ind_cols.size(); r++) { // loop over columns of graphe to get relevant row
                  if(ind_rows[r]==1) {
                    n_m1(rowc,_) = graphe(r,_);
                    rowc++;
                  }
                }
                
                NumericVector n_v1 = n_m1(0,_); //cut out level/state; (but there is only 1 row anyways ...)
                NumericVector n_v2;
                n_v2 = n_v1[ind_cols==1]; //cut out columns (node 2); n_v2 is now 1 parameter(number)
                
                NumericVector dummy_postcon; //apparently necessary to get from vector->double
                dummy_postcon = n_v2 * Data(p,node2); // interaction term of node2 (continuous variable); Data always contains the latest value
                potcon[con] = dummy_postcon[0]; // interaction term of node2 (continuous variable)
                con++;
                
                
              } // end if: node2: cat vs. cont
              
            } // end if: node2=node
            
          } // end: loop: other nodes (in node=continuous)
          
          double natpar;
          natpar = thresh_m(node,0) + sum(potcat) + sum(potcon);
          
          //limiting mean values to avoid infinite values in case the model is misspecified
          double bas;
          bas = 10;
          double po;
          po= 300;
          double epsi;
          epsi = pow(bas,po); 
          
          if(type_c[node]==2) { //gauss
            if (natpar>(epsi)) Rcpp::stop("Value of Gaussian node approaches Inf (> 10^300) in Gibbs sampler. Gaussian Submatrix (Covariance matrix) is not positive definite.");
            if (natpar<(-epsi)) Rcpp::stop("Value of Gaussian node approaches -Inf (< -10^300) in Gibbs sampler. Gaussian Submatrix (Covariance matrix) is not positive definite.");
            Data(p,node) = R::rnorm(natpar,1);
          } else if(type_c[node]==3) { //pois
            if (exp(natpar)<=0) Rcpp::stop("Lambda <= 0 for poisson node. Check the sign of the specified means and edge weights.");
            if (exp(natpar)>(epsi)) Rcpp::stop("Value of Poisson node approaches Inf (> 10^300) in Gibbs sampler. Check the sign of the specified means and edge weights.");
            Data(p,node) = R::rpois(exp(natpar));
          } else if(type_c[node]==4) { //exp
            if (natpar<=0) Rcpp::stop("Lambda(scale) <= 0 for exponential node. Check the sign of the specified means and edge weights.");
            if (natpar>(epsi)) Rcpp::stop("Value of Exponential node approaches Inf (> 10^300) in Gibbs sampler.  Check the sign of the specified means and edge weights.");
            Data(p,node) = R::rexp(1/natpar); // Rccp parameterization of rexp: scale instead of rate
          }
          
          
        } // end: if categorical vs. continuous
        
      } // end: loop nodes
      
    } // end: loop interactions
    
  } // end: loop cases
  
  
  return Data;
  
} // end of function


