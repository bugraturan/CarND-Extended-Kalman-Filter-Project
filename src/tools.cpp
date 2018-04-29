#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
	* Calculate the RMSE here.
  */
 	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
	* Calculate a Jacobian here.
	*/

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
			-(py/c1), (px/c1), 0, 0,
			py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}

VectorXd Tools::MapCartToPolar(const VectorXd& x_state) {
	/*
	This method maps the cart coordinates px, py, vx, and vy th rho, phi, and rho_dot
	*/

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float rho=0;
	float phi=0;
	float rho_dot=0;

	rho=sqrt(px*px+py*py);
	

	if (fabs(rho) < 0.0001) {
		phi = 0;
    	rho_dot = 0;
  	} else {
		phi=std::atan2(py, px);
    	rho_dot=(px*vx+py*vy)/rho;
	}
	
	VectorXd out = VectorXd(3);
	out(0)=rho;
	out(1)=phi;
	out(2)=rho_dot;

	return out;
 }

VectorXd Tools::MapPolarToCart(const VectorXd& x_state) {
	/*
	This method maps the polar coordinates rho, phi, and rho_dot to px, py, vx=0, and vy=0
	*/
	float px = 0;
	float py = 0;
	float vx = 0;
	float vy = 0;

	float rho=x_state(0);
	float phi=x_state(1);

	px = rho * cos(phi);
    py = rho * sin(phi);

	VectorXd out = VectorXd(1, 4);

	out(0)=px;
	out(1)=py;
	out(2)=vx;
	out(3)=vy;

	return out;
 }

double Tools::NormalizeAngle(double angle) {
	/*
	This method normalizes angles to be between -pi and pi
	*/
    double a = fmod(angle + M_PI, 2 * M_PI);
    return a >= 0 ? (a - M_PI) : (a + M_PI);
}