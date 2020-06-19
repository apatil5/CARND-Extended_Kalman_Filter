#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
 
  for(int i=0; i < estimations.size(); ++i){

    VectorXd diff = estimations[i] - ground_truth[i];

    //dot product
    VectorXd  diff_sq= diff.array() * diff.array();

    rmse +=diff_sq;
  }

  //Mean
  rmse = rmse / estimations.size();

  //RMS
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  float rho2 = x_state(0)*x_state(0) + x_state(1)*x_state(1);
  float rho  = sqrt(rho2);
  float rho3 = rho2*rho;
  float a11    = x_state(0)/rho;
  float b11    = x_state(1)/rho;
  float b21    = -(x_state(1)/rho2);
  float b22    = x_state(0)/rho2;
  float c      = (x_state(2)*x_state(1) - x_state(3)*x_state(0))/rho3;
  float c31    = x_state(1)*c;
  float c32    = -x_state(0)*c;
  float tv     = 1;
  //check ision by zero
  if(fabs(b21) < 0.0009 ){
    //rho =EPS;
    tv=100000;
  }
 
  //compute the Jacobian matrix
  Hj <<   a11,      b11,   0,     0,
          b21,   b22,   0,     0,
          c31,   c32,     a11,   b11;

  return tv*Hj;
}
