#include "kalman_filter.h"
#include "math.h"
#include<iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; 
  P_ = P_in; 
  F_ = F_in; 
  H_ = H_in; 
  R_ = R_in; 
  Q_ = Q_in; 
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  MatrixXd D   = H_ * P_ * H_.transpose() + R_;
  MatrixXd K   = P_*H_.transpose() * D.inverse();
  
  x_ = x_ + (K * (z-H_*x_));

  MatrixXd I = MatrixXd::Identity(4,4);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double rho   = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  
  double phi   = atan2(x_(1), x_(0));
  while (phi < -M_PI){phi +=M_PI*2;}
  while (phi>M_PI){phi-=2*M_PI;}
   
  double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
  VectorXd hx = VectorXd(3);
  hx << rho, phi, rho_dot;

  MatrixXd D   = H_ * P_ * H_.transpose() + R_;
  MatrixXd K   = P_*H_.transpose() * D.inverse();
  
  x_ = x_ + (K * (z-hx));

  MatrixXd I = MatrixXd::Identity(4,4);
  P_ = (I - K * H_) * P_;

}

