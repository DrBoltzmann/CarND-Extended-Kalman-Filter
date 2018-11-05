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
  
  // predict the state
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  
  std::cout<<"Predicted X is "<<x_<<std::endl;
  std::cout<<"Predicted P is "<<P_<<std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  // update the state by using Kalman Filter equations
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
  //int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // * update the state by using Extended Kalman Filter equations
  
  VectorXd h = VectorXd(3);
  
  h(0) = sqrt( x_(0) *  x_(0) + x_(1) * x_(1) );
  h(1) = atan2( x_(1) , x_(0) );
  h(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / h(0);  
  
  /*
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  */
  // Equations for h_func below
  //float eq1 = sqrt(px * px + py * py);
  
  //check division by zero
  /*
  if(eq1 < .00001) {
    px += .001;
    py += .001;
    eq1 = sqrt(px * px + py * py);
  }
  float eq2 = atan2(py,px);
  float eq3 = (px*vx+py*vy)/eq1;
  */
  //Feed in equations above
  //VectorXd H_func(3);
  
  //H_func << eq1, eq2, eq3;
    
  VectorXd y = z - h;

  y(1) = atan2(sin(y(1)),cos(y(1)));
/*  
  // Normalize the angle
  while (y(1)>M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1)<-M_PI) {
    y(1) += 2 * M_PI;
  }
*/
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  std::cout<<"Kalman Gain is " << K<<std::endl;
  
  // new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
