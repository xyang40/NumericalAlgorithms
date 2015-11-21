#include <Rcpp.h>
#include <stdint.h>
using namespace Rcpp;
using namespace std;
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define SQ(x) ((x)*(x))


// [[Rcpp::export]]
double aux1(NumericMatrix X,NumericVector cur_beta,int col,double rho,NumericMatrix D){
  int nrow = X.nrow();
  int ncol = X.ncol();
  
  double sum1=0;
  double sum3=0;
  NumericVector XBeta(ncol);
  NumericVector DBeta(ncol);                          
  for(register int i=0;i<ncol;i++){
    for(register int j=0;j<ncol;j++){
      XBeta[i] += X(i,j)*cur_beta[j];
      DBeta[i] += D(i,j)*cur_beta[j];
    }
    sum1 += SQ(XBeta[i])/(2*nrow);
    sum3 += rho*cur_beta[i]*DBeta[i]/2;
  }
  
  return sum1 - cur_beta[col] + sum3;
}

// [[Rcpp::export]]
double getVectorAbsMax(NumericVector nv){
   double max = 0.0;
   for(register int i=0;i<nv.size();i++){
     max = MAX(abs(nv[i]),max);
   }  
   return max;
}

// [[Rcpp::export]]
double getVectorAbsMax2(NumericVector nv,NumericVector aux){
   double max = 0.0;
   for(register int i=0;i<nv.size();i++){
     max = MAX(aux[i]*abs(nv[i]), max);
   }  
   return max;
}

// [[Rcpp::export]]
NumericMatrix bigscio_cd(const NumericMatrix X,NumericMatrix beta,const int maxit,const double stopcrit,
                         const double rho,const NumericMatrix D, const NumericVector lambda, 
                         const double nu,const double sigma,const double tau){
                           
                          const int nrow = X.nrow();
                          const int ncol = X.ncol();
                          
                          NumericVector Sig_diag(ncol);
                          for(register int i=0;i<ncol;i++){
                            for(register int j=0;j<nrow;j++){
                              Sig_diag[i] += SQ(X(j,i))/nrow;
                            }
                            Sig_diag[i] += rho*D(i,i);
                          }

                          for(register int col=0;col<ncol;col++){
                            
                            register int iter = 0;
                            while(iter<maxit){
                              
                              NumericVector XBeta(nrow);
                              for(register int i=0;i<nrow;i++){
                                for(register int j=0;j<ncol;j++){
                                  XBeta[i] += X(i,j)*beta(j,col);
                                }
                              }
                              NumericVector SigBeta(ncol);
                              NumericVector d(ncol);
                              for(register int i=0;i<ncol;i++){
                                for(register int j=0;j<nrow;j++){
                                  SigBeta[i] += XBeta[j]*X(j,i)/nrow; 
                                }
                                for(register int k=0;k<ncol;k++){
                                  SigBeta[i] += rho*D(i,k)*beta(k,col);
                                }
                                d[i] = (i==col?SigBeta[i]-1:SigBeta[i]);                         
                              }

                              NumericVector Sigd(ncol);
                              NumericVector mu(ncol);
                              
                              register int iter_in = 0;
                              while(iter_in<maxit){

                                NumericVector Xd(nrow);
                                for(register int i=0;i<nrow;i++){
                                  for(register int j=0;j<ncol;j++){
                                    Xd[i] += X(i,j)*d[j];
                                  }
                                }
                                
                                NumericVector Sigd(ncol);
                                double max_abs_mu = 0;
                                for(register int i=0;i<ncol;i++){
                                  for(register int j=0;j<nrow;j++){
                                    Sigd[i] += Xd[j]*X(j,i)/nrow; 
                                  }
                                  for(register int k=0;k<ncol;k++){
                                    Sigd[i] += rho*D(i,k)*d[k];
                                  }
   
                                  double a = Sig_diag[i];
                                  double b = (i==col?SigBeta[i]-1+Sigd[i]:SigBeta[i]+Sigd[i]);
                                  double c = beta(i,col) + d[i];
                                  if(c-(b/a)>0){
                                    mu[i] = -c + MAX(abs(c-(b/a)) - lambda[i]/a,0);
                                  }
                                  else if(c-(b/a)==0){
                                    mu[i] = -c;
                                  }
                                  else{
                                    mu[i] = -c - MAX(abs(c-(b/a)) - lambda[i]/a,0);
                                  }
                                  max_abs_mu = MAX(abs(mu[i]),max_abs_mu);
                                }
                                                 
                                if(max_abs_mu < stopcrit){
                                  break;
                                } 
                                else{
                                  d = d + mu;
                                }  
                               
                                iter_in++;
                             }
                             
                             double sum1 = 0;
                             double sum2 = 0;
                             for(register int i=0;i<ncol;i++){
                               sum1 += SQ(d[i]);
                               sum2 += tau*Sigd[i]*d[i];
                             }
                             double alpha = 1;
                             double delta = sum1+sum2+getVectorAbsMax2(beta(_,col)+d,lambda)+getVectorAbsMax2(beta(_,col),lambda);
                             double t1 = aux1(X, beta(_,col)+alpha*d,col,rho,D) + getVectorAbsMax2(beta(_,col)+alpha*d,lambda);
                             double t2 = aux1(X, beta(_,col),col,rho,D) + getVectorAbsMax2(beta(_,col),lambda) + alpha*sigma*delta;
                            
                             register int k = 0;
                             while(t1 > t2){
                               k++;
                               alpha = pow(nu,k);
                               t1 = aux1(X, beta(_,col)+alpha*d,col,rho,D) + getVectorAbsMax2(beta(_,col)+alpha*d,lambda);
                               t2 = aux1(X, beta(_,col),col,rho,D) + getVectorAbsMax2(beta(_,col),lambda) + alpha*sigma*delta;
                             }
                             
                             if(getVectorAbsMax(alpha*d) < stopcrit){
                                break;
                             } 
                             else{
                              for(register int i=0;i<ncol;i++){
                                if((lambda[i]==0)||(abs(beta(i,col))>stopcrit||(abs(d[i]))>lambda[i])){    
                                  beta(i,col) += alpha*d[i];
                                }
                              } 
                             }                
                             
                             iter++;
                            }                          
                          }
                          
                          for(register int i=0;i<ncol;i++){
                            for(register int j=0;j<ncol;j++){
                              if(abs(beta(i,j))<=abs(beta(j,i))){
                                beta(j,i)=beta(i,j);
                              }
                              else{ 
                                beta(i,j)=beta(j,i);
                              }
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
lambda<-rep(1,p)
maxit=1000
stopcrit=10^-6

nu=0.5
sigma=0.5
tau=0

solve(cov(X))

ptm <- proc.time()
res_xi_cd=bigscio_cd(X,beta=matrix(0.1,p,p)+diag(0.5,p),maxit,stopcrit,rho,D,lambda=rep(0.01,p), nu,sigma,tau)
proc.time() - ptm   

*/
