#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define sq(x) ((x)*(x))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// [[Rcpp::export]]
NumericMatrix bigscio_cg(const NumericMatrix X,NumericMatrix beta,const int maxit,const double stopcrit,
                         const double rho,const NumericMatrix D, const NumericVector lambda){
                           
                           const int nrow = X.nrow();
                           const int ncol = X.ncol();
 
                           for(register int col=0;col<ncol;col++){
                             
                              NumericVector XBeta(nrow);
                              for(register int i=0;i<nrow;i++){
                                for(register int j=0;j<ncol;j++){
                                  XBeta[i] += X(i,j)*beta(j,col);
                                }
                              }
                              NumericVector SigBeta(ncol);
                              NumericVector r(ncol);
                              NumericVector p(ncol);
                              for(register int i=0;i<ncol;i++){
                                for(register int j=0;j<nrow;j++){
                                  SigBeta[i] += XBeta[j]*X(j,i)/nrow; 
                                }
                                for(register int k=0;k<ncol;k++){
                                  SigBeta[i] += rho*D(i,k)*beta(k,col);
                                }                  
                                r[i] = (i==col?1-SigBeta[col]:-SigBeta[i]);
                                p[i] = r[i];
                              }
                   
                              register int iter = 0;
                              while(iter<maxit){                         
                                NumericVector XBeta(nrow);
                                for(register int i=0;i<nrow;i++){
                                  for(register int j=0;j<ncol;j++){
                                    XBeta[i] += X(i,j)*beta(j,col);
                                  }
                                }
                                NumericVector SigBeta(ncol);
                                for(register int i=0;i<ncol;i++){
                                  for(register int j=0;j<nrow;j++){
                                    SigBeta[i] += XBeta[j]*X(j,i)/nrow; 
                                  }
                                  for(register int k=0;k<ncol;k++){
                                    SigBeta[i] += rho*D(i,k)*beta(k,col);
                                  }                  
                                }

                                double left_sum = 0;
                                NumericVector Xp(nrow);
                                NumericVector XXp(ncol);
                                for(register int i=0;i<nrow;i++){
                                  for(register int j=0;j<ncol;j++){
                                    Xp[i] += X(i,j)*p[j];
                                  }
                                  left_sum += sq(Xp[i])/nrow;
                                }
                                double right_sum = 0;
                                double numerator=0.0;
                                NumericVector Dp(ncol);
                                for(register int i=0;i<ncol;i++){
                                  for(register int j=0;j<nrow;j++){
                                    XXp[i] += Xp[j]*X(j,i)/nrow;
                                  }
                                  numerator+=sq(r[i]);
                                  for(register int j=0;j<ncol;j++){
                                    Dp[i] += D(i,j)*p[j];
                                  }
                                  right_sum += rho*p[i]*Dp[i];
                                }
                                NumericVector SigP = XXp + rho*Dp;
                        
                                double alpha = numerator/(left_sum+right_sum);                             
                                for(register int i=0;i<ncol;i++){
                                 if(lambda[i]==0 || abs(beta(i,col) > stopcrit || abs((i==col?abs(SigBeta[i]-1):abs(SigBeta[i])))>lambda[i])){
                                     beta(i,col) += alpha*p[i];
                                 }
                                }                            
                               
                                NumericVector r_new = r - alpha*SigP;
                                double max_r = 0.0;
                                for(register int i=0;i<ncol;i++){
                                  max_r = MAX(abs(r_new[i]),max_r);
                                }                   
                                if(max_r < stopcrit){
                                  break;
                                }      
                                
                                double sum1=0.0;
                                double sum2=0.0;
                                for(register int i=0;i<ncol;i++){
                                  sum1+=sq(r_new[i]);
                                  sum2+=sq(r[i]);
                                }                            
                                p = r_new + (sum1/sum2)*p;
                                r = r_new;
                                iter++;
                              }
                            }
                          
                           for(register int i=0;i<ncol;i++){
                             for(register int j=0;j<ncol;j++){
                               if(abs(beta(i,j))<=abs(beta(j,i))){beta(j,i)=beta(i,j);}
                               else{ beta(i,j)=beta(j,i);}
                             }
                           }                     
                          return beta;
                         }
                         
                         
                         
                         
                         

   

/*** R

n<-100
p<-10
Omega<-matrix(0.5,p,p)+diag(0.5,p)
sv<-svd(Omega)
set.seed(20)
X<-matrix(rnorm(n*p),n,p)%*%(sv$u%*%diag(1/sqrt(sv$d))%*%t(sv$v))

rho<-0
D<-diag(rep(1,p))
lambda<-rep(0,p)
maxit = 1000
solve(cov(X))
stopcrit = 10^-6

ptm <- proc.time()
res_xi_cg=bigscio_cg(X,beta=matrix(0,p,p)+diag(0.5,p),maxit,stopcrit,rho,D,lambda)
proc.time() - ptm
*/



