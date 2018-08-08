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

  // initializing matrices
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
  
  // H is the measurement matrix that projects from the 4D state, your belief about the object's current state, into the 2D measurement space of the sensor (Lesson 5 - 11. Laser Measurements)
  H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;
  
  // * Finish initializing the FusionEKF
  // * Set the process and measurement noises

  ekf_.F_ = MatrixXd(4, 4); // 4 x 4 state transition matrix
  ekf_.F_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  ekf_.P_ = MatrixXd(4, 4); // 4 x 4 state covariance matrix
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  // acceleration noise components
   noise_ax = 9;
   noise_ay = 9;

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
    // * Initialize the state ekf_.x_ with the first measurement.
    // * Create the covariance matrix.
    // * Remember: you'll need to convert radar from polar to cartesian coordinates.
    
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar (theta, r) to cartesian (x, y) coordinates and initialize state.
      ekf_.x_(0) = ro*cos(theta);
      ekf_.x_(1) = ro*sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize the state
      ekf_.x_(0) << x;
      ekf_.x_(1) << y;
    }

    // done initializing, no need to predict or update
    // ekf_.F << set to 1 diagonal matrix
    //previous timestamp = measurement_pack.timestamp;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */ 
  
  // define dt in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  // Modify the F matrix to integrate the time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // Set the process covariance matrix Z
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay, 0,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
    
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // * Use the sensor type to perform the update step.
  // * Update the state and covariance matrices.
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); //call the function defined in tools.cpp to calculate the Jacobian
    ekf_.R_ = R_radar_;
    ekf_.UdateEKF(measurement_pack.raw_measurements_);
  
  } else {
   // Laser updates
  ekf_H._ = H_laser_;
  ekf_R._ = R_laser_;
  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
