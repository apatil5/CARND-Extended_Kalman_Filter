#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
 
  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;
    
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      double theta = measurement_pack.raw_measurements_[1];
      ekf_.x_ << measurement_pack.raw_measurements_[0]*cos(theta), 
                 measurement_pack.raw_measurements_[0]*sin(theta), 
                 0, 
                 0; 
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
                 measurement_pack.raw_measurements_[1], 
                 0, 
                 0; 
    }
    
    previous_timestamp_=measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
 }

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  double dt2=dt*dt;
  double dt3=dt2*dt;
  double dt4=dt3*dt;
  
  double noise_ax_sq=81;
  double noise_ay_sq=81;
  
  
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  float q33 = noise_ax_sq*dt2;
  float q44 = noise_ay_sq*dt2;
  float q11= q33*dt2/4;
  float q13= q33*dt/2;
  float q22= q44*dt2/4;
  float q24= q44*dt/2;

  ekf_.Q_ << q11, 0,   q13,  0,
             0,   q22,   0, q24,
             q13,  0,   q33,  0,
             0,   q24,   0, q44;

  ekf_.Predict();

  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      	
     Hj_= tools.CalculateJacobian(ekf_.x_);
     ekf_.Init( ekf_.x_,ekf_.P_,ekf_.F_,Hj_,R_radar_,ekf_.Q_);
     ekf_.UpdateEKF(measurement_pack.raw_measurements_);   

  } else {
  
     ekf_.Init(  ekf_.x_,ekf_.P_,ekf_.F_,H_laser_,R_laser_,ekf_.Q_);
     ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout<<"H"<<ekf_.H_<<endl;
}
