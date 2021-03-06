#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method that calculates h(x').
  */
  VectorXd MapCartToPolar(const VectorXd& x_state);

  /**
  * A helper method that calculates inverse of h(x').
  */
  VectorXd MapPolarToCart(const VectorXd& x_state);

  /**
  * A helper method that normalizes angles between -pi and pi.
  */
  double NormalizeAngle(double angle);

};

#endif /* TOOLS_H_ */
