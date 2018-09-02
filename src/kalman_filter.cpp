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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  // update the state by using Kalman Filter equations
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  //MatrixXd PHt = P_ * Ht;
  MatrixXd K =  P_ * Ht * Si;
  
  // New state
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
  //KF Prediction step
  //x = F * x + u;
  //MatrixXd Ft = F.transpose();
  //P = F * P * Ft + Q;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // * update the state by using Extended Kalman Filter equations
  
  float x = x_(0);
  float y = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float ro = sqrt(x * x + y * y);
  float theta = atan2(y, x);
  //float rho_dot = (x * vx + y * vy) / rho;
  
  float ro_dot;  
  if (fabs(ro) < 0.0001) {
    ro_dot = 0;
  } else {
    ro_dot = (x*vx + y*vy)/ro;
  }
  
  //VectorXd z_pred = VectorXd(3);
  VectorXd z_pred(3);
  
  z_pred << ro, theta, ro_dot;
  
  VectorXd yy = z - z_pred;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = P_ * Ht * Si;
  
  // new state
  x_ = x_ + (K * yy);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
