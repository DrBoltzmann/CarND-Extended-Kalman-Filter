#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  // predict the state
  //x = F * x + u; where u is external motion
  x = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // update the state by using Kalman Filter equations
  VectorXd y = z - H * x;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P * Ht * Si;
  
  // New state
  x = x + (K * y);
  P = (I - K * H) * P;
  
  //KF Prediction step
  //x = F * x + u;
  //MatrixXd Ft = F.transpose();
  //P = F * P * Ft + Q;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // * update the state by using Extended Kalman Filter equations
  
  float x = ekf_.x_(0);
  float y = ekf_.x_(1);
  float vx = ekf_.x_(2);
  float vy = ekf_.x_(3);
  
  float rho = sqrt(x * x + y * y);
  float theta = atan2(y, x);
  float ro_dot = (x * vx + y * vy) / rho;
  
   VectorXd z_pred = VectorXd(3);
   z_pred << rho, theta, ro_dot;
  
  VectorXd y = z - z_pred;
  
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = P * Ht * Si;
  
  // new state
  
  x = x + (K * y);
  P = (I - K * H) * P; 
  
  // KF Prediction step
  
  x = F * x + u;
  MatrixXd Ft = F.transpose();
  P = F * P * Ft + Q;
  
}
