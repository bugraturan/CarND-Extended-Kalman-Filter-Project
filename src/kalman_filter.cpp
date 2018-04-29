#include <iostream>
#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  /* 
  The prediction step is the same for the KF and EKF
  */
 	x_          = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_          = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /* 
  For LASER measurements we can implement the same process as in the lecture
  */

 	VectorXd z_pred = H_ * x_;
	VectorXd y      = z - z_pred;
	MatrixXd Ht     = H_.transpose();
	MatrixXd S      = H_ * P_ * Ht + R_;
	MatrixXd Si     = S.inverse();
	MatrixXd PHt    = P_ * Ht;
	MatrixXd K      = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /* 
    For RADAR measurements we need to use a different update process
  */
  Tools tools; // Helper methods

  VectorXd z_pred = tools.MapCartToPolar(x_); // Convert from cart to polar

  VectorXd y      = z - z_pred;
  y(1) = tools.NormalizeAngle(y(1)); // We need to be sure that the resulting phi is between -pi and pi

  // The rest is similar to the KF
  MatrixXd Ht     = H_.transpose();
  MatrixXd S      = H_ * P_ * Ht + R_;
  MatrixXd Si     = S.inverse();
  MatrixXd PHt    = P_ * Ht;
  MatrixXd K      = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
