#include <TMB.hpp> // Links in the TMB libraries

#include "sim_multinom.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;
  
  DATA_MATRIX(H); // The halibut data for 300 hooks
  DATA_MATRIX(A); // The halibut data for 700 hooks
  DATA_MATRIX(X); // Covariates matrix 
  PARAMETER_VECTOR(betat); // Parameters for lambda t (target species)
  PARAMETER_VECTOR(betant); // Parameters for lambda nt (non-target species)
  PARAMETER(theta); // Parameters for pnt (probability of caught non-target fish escaping)
  DATA_VECTOR(s); // Soak time

  int n = H.col(0).size(); // Length of the column of matrix H
  int n_k = H.row(0).size(); // Length of the row of matrix H (4)
  int n2 = A.col(0).size(); // Length of the column of matrix A
  int n_k2 = A.row(0).size(); // Length of the row of matrix A (3)
  PARAMETER_VECTOR(omegat); // The random field for lambda t
  PARAMETER_VECTOR(omegant); // The random field for lambda nt
  
  DATA_MATRIX(D); // Distance matrix
  PARAMETER(lognut);
  PARAMETER(lognunt);
  PARAMETER(logPhit);
  PARAMETER(logPhint);
  PARAMETER(logSigmat);
  PARAMETER(logSigmant);
  // spatial parameters
  Type nut=exp(lognut); // smoothness parameter 
  Type nunt=exp(lognunt);
  Type phit=exp(logPhit); // range parameter 
  Type phint=exp(logPhint); 
  Type sigt=exp(logSigmat); // variance
  Type signt=exp(logSigmant);
  
  vector <Type> nll_joint(4); nll_joint.setZero();
  
  // SIMULATE {
  //   //SIMULATE A DISTANCE MATRIX
  //   for (int i = 0; i < n; i++) {
  //     for (int j = 0; j < n; j++) {
  //       if (i == j) D(i,j) = Type(0.0);
  //       if (i < j) D(i,j) = exp(rnorm(Type(0.28),Type(0.3)));
  //       if (i > j) D(i,j) = D(j,i);
  //     }
  //   }
  //   REPORT(D);
  // }
  
  // Covaraince matirx for random field for lambda t
  matrix<Type> St(n,n);
  St.setZero();
  for(int i=0; i<n; ++i){
    St(i,i) = sigt*sigt;
  }
  for(int i=0; i<n; ++i){
    for(int j=i+1; j<n; ++j){
      St(i,j) = sigt*sigt*matern(D(i,j), phit, nut); //Matern Covariance Function
      St(j,i) = St(i,j);
    }
  }
  
  // Covaraince matirx for random field for lambda nt
  matrix<Type> Snt(n,n);
  Snt.setZero();
  for(int i=0; i<n; ++i){
    Snt(i,i) = signt*signt;
  }
  for(int i=0; i<n; ++i){
    for(int j=i+1; j<n; ++j){
      Snt(i,j) = signt*signt*matern(D(i,j), phint, nunt); //Matern Covariance Function
      Snt(j,i) = Snt(i,j);
    }
  }
  
  REPORT(St);
  REPORT(Snt);
  
  nll_joint(0) += density::MVNORM(St)(omegat);
  nll_joint(1) += density::MVNORM(Snt)(omegant);
    
    SIMULATE{
      // matrix<Type> St(n,n);
      // St.setZero();
      for(int i=0; i<n; ++i){
        St(i,i) = sigt*sigt;
      }
      for(int i=0; i<n; ++i){
        for(int j=i+1; j<n; ++j){
          St(i,j) = sigt*sigt*matern(D(i,j), phit, nut); //Matern Covariance Function
          St(j,i) = St(i,j);
        }
      }
      
      // Covaraince matirx for random field for lambda nt
      // matrix<Type> Snt(n,n);
      // Snt.setZero();
      for(int i=0; i<n; ++i){
        Snt(i,i) = signt*signt;
      }
      for(int i=0; i<n; ++i){
        for(int j=i+1; j<n; ++j){
          Snt(i,j) = signt*signt*matern(D(i,j), phint, nunt); //Matern Covariance Function
          Snt(j,i) = Snt(i,j);
        }
      }
      REPORT(St);
      REPORT(Snt);
      
      density::MVNORM_t <Type> t_mvnorm(St);
      density::MVNORM_t <Type> nt_mvnorm(Snt);
      
      vector <Type> temp_omegat(n);
      vector <Type> temp_omegant(n);
      
      t_mvnorm.simulate(temp_omegat);
      nt_mvnorm.simulate(temp_omegant);
      
      omegat = temp_omegat;
      omegant = temp_omegant;
      
      REPORT(omegat);
      REPORT(omegant);
    
      lognut = log(nut);
      lognunt = log(nunt);
      logPhit = log(phit);
      logPhint = log(phint);
      logSigmat = log(sigt);
      logSigmant = log(signt);
      REPORT(lognut);
      REPORT(lognunt);
      REPORT(logPhit);
      REPORT(logPhint);
      REPORT(logSigmat);
      REPORT(logSigmant);
    }
    
    SIMULATE{
      for (int i = 0; i < n; i++){
        s(i) = exp(rnorm(log(Type(450)),Type(0.2)));
      }
      REPORT(s);
    }
    
    // ldat - the relative abundance indices for target species (halibut)
    vector<Type> ldat(n);     
    for(int i=0; i<n; ++i){
      // Adding random field and covariates to lambda t
      ldat(i)=exp((vector<Type>(X.row(i))*betat).sum()+omegat(i));
    }
    ADREPORT(ldat);
    
    // ldant - the relative abundance indices for non-target species 
    vector<Type> ldant(n);     
    for(int i=0; i<n; ++i){
      // Adding omega (random field) and covariates to lambda nt 
      ldant(i)=exp((vector<Type>(X.row(i))*betant).sum()+omegant(i));
    }
    ADREPORT(ldant);
    
    // pnt - the probability of caught non-target species escaping 
    vector<Type> pnt(n);     
    for(int i=0; i<n; ++i){
      pnt(i)=invlogit(theta);
    }
    ADREPORT(pnt);
    
    // Corresponding probabilities for N_b, N_t, N_nt and N_e for 300 hooks
    // N_b is the nmuber of baited hooks
    // N_t is the number of target species caught
    // N_nt is the number of non-target species caught
    // N_e is the number of empty hooks
  matrix<Type> p(n, n_k); 
  for(int i=0;i<n;i++){
    p(i,0) = exp(-(ldat(i)+ldant(i))*s(i));
    p(i,1)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldat(i)/(ldat(i)+ldant(i)));
    p(i,2)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldant(i)/(ldat(i)+ldant(i)))*(Type(1)-pnt(i));
    p(i,3)=1-p(i,0)-p(i,1)-p(i,2);
  }
  
  SIMULATE{
    for (int i=0; i<n; i++){
      X(i,0) = Type(1.0);
    }
    REPORT(X);
    
    // vector<Type> ldat(n);     
    for(int i=0; i<n; ++i){
      // Adding random field and covariates to lambda t
      ldat(i)=exp((vector<Type>(X.row(i))*betat).sum()+omegat(i));
    }
    REPORT(ldat);
    REPORT(betat);
    
    // ldant - the relative abundance indices for non-target species 
    // vector<Type> ldant(n);     
    for(int i=0; i<n; ++i){
      // Adding omega (random field) and covariates to lambda nt 
      ldant(i)=exp((vector<Type>(X.row(i))*betant).sum()+omegant(i));
    }
    REPORT(ldant);
    REPORT(betant);
    
    for(int i=0; i<n; ++i){
      pnt(i)=invlogit(theta);
    }
    REPORT(pnt);
    REPORT(theta);
    
    for(int i=0;i<n;i++){
      p(i,0) = exp(-(ldat(i)+ldant(i))*s(i));
      p(i,1)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldat(i)/(ldat(i)+ldant(i)));
      p(i,2)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldant(i)/(ldat(i)+ldant(i)))*(Type(1)-pnt(i));
      p(i,3)=1-p(i,0)-p(i,1)-p(i,2);
    }
    REPORT(p);
  }
  
  // Corresponding probabilities for N_t, N_nt and N_b+e for 700 hooks
  // N_t is the number of target species caught
  // N_nt is the number of non-target species caught
  // N_b+e is the number of empty hooks (baited + unbaited)
  matrix<Type> q(n2, n_k2); 
  for(int i=0;i<n;i++){
    q(i,0)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldat(i)/(ldat(i)+ldant(i)));
    q(i,1)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldant(i)/(ldat(i)+ldant(i)))*(Type(1)-pnt(i));
    q(i,2)=1-q(i,0)-q(i,1);
  }
  
  SIMULATE{
    for(int i=0;i<n;i++){
      q(i,0)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldat(i)/(ldat(i)+ldant(i)));
      q(i,1)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldant(i)/(ldat(i)+ldant(i)))*(Type(1)-pnt(i));
      q(i,2)=1-q(i,0)-q(i,1);
    }
    REPORT(q);
  }
  
  // The likelihood function for multinomial distribution
  for(int i=0; i < n; i++){
    vector<Type> p_row=p.row(i);
    vector<Type> H_row=H.row(i);
    vector<Type> q_row=q.row(i);
    vector<Type> A_row=A.row(i);
    nll_joint(2) -= dmultinom(H_row,p_row,true);
    nll_joint(2) -= dmultinom(A_row,q_row,true);
  }
  
  SIMULATE{
    for (int i = 0; i < n; i++){
      vector <Type> p_row=p.row(i);
      vector<Type> H_row=H.row(i);
      vector<Type> q_row=q.row(i);
      vector<Type> A_row=A.row(i);
      
      H_row = sim_multinom(Type(300.0),p_row);
      A_row = sim_multinom(Type(700.0),q_row);
      
      H.row(i)=H_row;
      A.row(i)=A_row;
    }
    REPORT(p);
    REPORT(H);
    REPORT(q);
    REPORT(A);
  }
  
  REPORT(nll_joint);
  
  REPORT(ldat);
  REPORT(ldant);
  REPORT(pnt);
  REPORT(p);
  REPORT(betat);
  REPORT(betant);
  REPORT(omegat);
  REPORT(omegant);
  ADREPORT(phit);
  ADREPORT(phint);
  Type vart = sigt*sigt;
  Type varnt = signt*signt;
  ADREPORT(vart);
  ADREPORT(varnt);
  
  Type nll = nll_joint.sum();
  return nll;
}
