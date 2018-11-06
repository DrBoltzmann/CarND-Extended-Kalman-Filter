#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using std::vector;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; // object state
  P_ = P_in; // object covariance matrix
  F_ = F_in; // state transition matrix
  H_ = H_in; // measurement matrix
  R_ = R_in; // measurement covariance matrix
  Q_ = Q_in; // process covariance matrix
}

void KalmanFilter::Predict() {
  
  std::cout<<"Initial X is "<<x_<<std::endl;
  std::cout<<"Initial P is "<<P_<<std::endl;
  std::cout<<"Initial Q is "<<Q_<<std::endl;  
  
  // Predict the state
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  
  std::cout<<"Predicted X is "<<x_<<std::endl;
  std::cout<<"Predicted P is "<<P_<<std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  // Update the state by using Kalman Filter equations
  VectorXd z_pred = H_ * x_;
  MatrixXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K =  PHt * Si;
  
  // New state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // * update the state by using Extended Kalman Filter equations
  
  VectorXd h = VectorXd(3);
  
  h(0) = sqrt( x_(0) *  x_(0) + x_(1) * x_(1) );
  h(1) = atan2( x_(1) , x_(0) );
  h(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / h(0);  
    
  VectorXd y = z - h;

  y(1) = atan2(sin(y(1)),cos(y(1)));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  std::cout<<"Kalman Gain is " << K<<std::endl;
  
  // New state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
