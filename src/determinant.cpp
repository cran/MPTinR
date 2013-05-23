
// includes from the plugin
#include <RcppEigen.h>
#include <Rcpp.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP determinant( SEXP SS, SEXP Mcast, SEXP Acast, SEXP Bcast, SEXP Ccast, SEXP patternCast, SEXP ineq, SEXP samples) ;
}

// definition

SEXP determinant( SEXP SS, SEXP Mcast, SEXP Acast, SEXP Bcast, SEXP Ccast, SEXP patternCast, SEXP ineq, SEXP samples){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::ArrayXd;
const int S = Rcpp::as<int>(SS);
const int M = Rcpp::as<int>(Mcast);
const int s = Rcpp::as<int>(samples);
//Rcout << "Ineq:\n"<<"\n";
const Map<MatrixXd> Ineq(as<Map<MatrixXd> >(ineq));
const Map<MatrixXi> A(as<Map<MatrixXi> >(Acast));
const Map<MatrixXi> B(as<Map<MatrixXi> >(Bcast));
const Map<VectorXd> c(as<Map<VectorXd> >(Ccast));
const Map<MatrixXd> pattern(as<Map<MatrixXd> >(patternCast));
//Rcout << "pattern:\n"<< pattern <<"\n";

RNGScope scope;


// for sampling random data
bool found = false;
VectorXd thetaTMP;
RowVectorXd theta;
VectorXd IneqT;
bool check1;
bool check2;

// other vars
MatrixXd Theta;
VectorXd p;
MatrixXd V;
VectorXd pc;
MatrixXd delta0;
double I;

// for output
NumericVector out(3);

out[0] = 0;
out[1] = 0;
out[2] = 0;

for (int i = 0; i < s; i++) {
  found = false;
  while (!found) {
    check1 = true;
    check2 = true;
    thetaTMP = Rcpp::as<VectorXd>(rbeta(S, 0.5, 0.5));
    theta = thetaTMP.transpose();
  
    IneqT = (theta*Ineq.transpose());
  
    for(int i = 0; i< IneqT.size(); i++){
        if(IneqT(i) < 0)
            check1 = false;
    }
    for(int i = 0; i<S;i++){
        if(theta(i) == 1 || theta(i) == 0)
            check2 = false;
    }
    if (check1 && check2) {
      found = true;
    }
    out[0] ++;
  }
  
  Theta = VectorXd::Ones(M)*theta;
  // Rcout << "Theta:" << Theta << "\n";
  //Rcout << "IneqT \\< 0:" << (IneqT < 0) << "\n";
  //Rcout << ((theta*Ineq.transpose())<0) << "\n";
  //Rcout << "thetaTMP:\n"<<thetaTMP<<"\n";
  //Rcout << "theta:\n"<<theta<<"\n";
  
  MatrixXd TA_tmp(Theta.rows(), Theta.cols());
  MatrixXd TB_tmp(Theta.rows(), Theta.cols());
  
  for (int r = 0; r < Theta.rows(); r++) {
      for (int c = 0; c < Theta.cols(); c++) {
          //Rcout << "Theta: " << Theta(r, c) << "   A: " << A(r, c) << "\n";
          TA_tmp(r, c) = std::pow(Theta(r, c), A(r, c));
          TB_tmp(r, c) = std::pow((1 - Theta(r, c)), B(r, c));
      }
  }
  
  // VectorXd TA_tmp2 = TA_tmp.rowwise().prod();
  // VectorXd TB_tmp2 = TB_tmp.rowwise().prod();
  p = TA_tmp.rowwise().prod().array() * TB_tmp.rowwise().prod().array() * c.array();
  
  //Rcout << "p diagonal:\n" << p.asDiagonal() << "\n";
  
  // MatrixXd P_tmp = p.asDiagonal();
  // MatrixXd V = pattern * P_tmp;
  V = pattern * p.asDiagonal();
  //Rcout << "V:\n"<< V <<"\n";
  //MatrixXd out = Theta.array().pow(A.array());
  
  pc = V.rowwise().sum();
  VectorXd D(pc.size());
  for (int c = 0; c < D.size(); c++) {
      if (pc(c) == 0)
          D(c) = 0;
      else D(c) = 1 / pc(c);
  }
  
  // MatrixXd thetaD = theta.asDiagonal();
  //Rcout << "thetaD:\n"<<thetaD<<"\n";
  // MatrixXd ABA = (A-(A+B)).cast<double>();
  // Rcout << "ABA:\n"<<ABA<<"\n";
  delta0 = V*(A.cast<double>()-(A+B).cast<double>()*theta.asDiagonal());
  
  I = (delta0.transpose() * D.asDiagonal() * delta0 * (theta.array()*(1 - theta.array())).matrix().asDiagonal().inverse() * M_PI * M_PI).determinant();
  // double detI = I.determinant();
  
  out[1] += I;
  out[2] += std::sqrt(I);
  
}
return wrap(out);

END_RCPP
}



