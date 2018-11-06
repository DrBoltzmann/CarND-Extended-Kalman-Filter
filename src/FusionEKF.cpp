#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // Initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Define laser covariance matrix
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Define radar covariance matrix
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // Define matrix for Jacobian
  Hj_<< 1,1,0,0,
        1,1,0,0,
        1,1,1,1;
  
  // * Finish initializing the FusionEKF
  // * Set the process and measurement noises

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1;  
  
  // H is the measurement matrix that projects from the 4D state, the belief about the object's current state, into the 2D measurement space of the sensor (Lesson 5 - 11. Laser Measurements)

  H_laser_ << 1,0,0,0,
              0,1,0,0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    
    // * Remember: you'll need to convert radar from polar to cartesian coordinates.
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    ekf_.P_ = MatrixXd(4, 4); // 4 x 4 state covariance matrix
    ekf_.P_ << 1,0,0,0,
               0,1,0,0,
               0,0,1,0,
               0,0,0,1;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize the state
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1); 
      ekf_.x_(2) = 0.0;
      ekf_.x_(3) = 0.0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar (theta, r) to cartesian (x, y) coordinates and initialize state.
      
      double ro = measurement_pack.raw_measurements_(0);
      double phi= measurement_pack.raw_measurements_(1);
      double ro_dot = measurement_pack.raw_measurements_(2);
      
      ekf_.x_(0) = ro * cos(phi);
      ekf_.x_(1) = ro * sin(phi);
      ekf_.x_(2) = ro_dot * cos(phi);
      ekf_.x_(3) = ro_dot * sin(phi);
    }
    
     ekf_.x_(2) = 5.2;
     ekf_.x_(3) = 1.8/1000.0;
    
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // * Update the state transition matrix F according to the new elapsed time.
  //    - Time is measured in seconds.
  //  * Update the process noise covariance matrix.
  //  * Use noise_ax = 9 and noise_ay = 9 for your Q matrix. 
  
  // Convert time to proper unit
  double dt = (measurement_pack.timestamp_ - previous_timestamp_ ) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  float noise_ax = 9.0;
  float noise_ay = 9.0; 

  double dt_2 = dt*dt;
  double dt_3 = dt_2*dt;
  double dt_4 = dt_2*dt_2;
  
  // Set the process covariance matrix Z
  ekf_.Q_ = MatrixXd(4, 4);
  
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;  
  
  // Modify the F matrix to integrate time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // * Use the sensor type to perform the update step.
  // * Update the state and covariance matrices.
  
  // Radar updates
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Hj_ = tools.CalculateJacobian(ekf_.x_); //call the function defined in tools.cpp to calculate the Jacobian
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }  else if
  
  // Laser updates
  (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
  ekf_.H_ = H_laser_;
  ekf_.R_ = R_laser_;
  ekf_.Update(measurement_pack.raw_measurements_);
  }
  
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
