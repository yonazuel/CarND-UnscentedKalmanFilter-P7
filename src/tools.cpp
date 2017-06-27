#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse = VectorXd(4);
	rmse << 0,0,0,0;

	for (int i = 0; i<estimations.size(); ++i){
		VectorXd res = estimations[i] - ground_truth[i];
		VectorXd se = res.array() * res.array();
		rmse += se;
	}

	rmse /= estimations.size();
	return rmse.array().sqrt();
}