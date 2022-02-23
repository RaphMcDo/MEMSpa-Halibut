#include <TMB.hpp> // Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x); // The halibut data 
  PARAMETER_VECTOR(theta); // Parameters
  DATA_VECTOR(s); // Soak time
  int n = x.col(0).size(); // Length of the column of matrix x 
  int n_k = x.row(0).size(); // Length of the row of matrix x (4)
  
  // ldat - the relative abundance index for target species (halibut)
  vector<Type> ldat(n);     
  for(int i=0; i<n; ++i){
    ldat(i)=exp(theta(0));
  }
  ADREPORT(ldat);
  
  // ldant - the relative abundance index for non-target species 
  vector<Type> ldant(n);     
  for(int i=0; i<n; ++i){
    ldant(i)=exp(theta(1));
  }
  ADREPORT(ldant);
  
  // pnt - the escaping rate for non-target species 
  vector<Type> pnt(n);     
  for(int i=0; i<n; ++i){
    pnt(i)=invlogit(theta(2));
  }
  ADREPORT(pnt);
  
  Type nll;
  
  matrix<Type> p(n, n_k); 
  for(int i=0;i<n;i++){
    p(i,0) = exp(-(ldat(i)+ldant(i))*s(i));
    p(i,1)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldat(i)/(ldat(i)+ldant(i)));
    p(i,2)=(Type(1)-exp(-(ldat(i)+ldant(i))*s(i)))*(ldant(i)/(ldat(i)+ldant(i)))*(Type(1)-pnt(i));
    p(i,3)=1-p(i,0)-p(i,1)-p(i,2);
  }
  
  for(int i=0; i < n; i++){
    vector<Type> p_row=p.row(i);
    vector<Type> x_row=x.row(i);
    nll -= dmultinom(x_row,p_row,true);
  }
  
  REPORT(ldat);
  REPORT(ldant);
  REPORT(pnt);
  return nll;
}
