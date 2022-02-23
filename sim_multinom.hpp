
using namespace density;
using namespace Eigen;

template<class Type>
vector <Type> sim_multinom(Type N, vector <Type> p) {
  int K = p.size();
  vector <Type> response(K); response.setZero();
  vector <Type> p_denom(K);
  vector <Type> min_N(K);
  // if(isDouble<Type>::value && of->do_simulate){
    for (int k = 0; k < K; k++) {
      if (k == 0) {p_denom(k) = 1; min_N(k) = 0;}
      if (k == 1) {p_denom(k) = 1- p(k-1); min_N(k) = response(k-1);}
      if (k == 2) {p_denom(k) = 1- p(k-2)-p(k-1); min_N(k) = response(k-2)+response(k-1);}
      if (k == 3) {p_denom(k) = 1- p(k-3)-p(k-2)-p(k-1); min_N(k) = response(k-3)+response(k-2)+response(k-1);}
      response(k) = rbinom(N-min_N(k), p(k)/p_denom(k));
    }
  // }
  return response;
}