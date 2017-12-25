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

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	VectorXd residual;

	for (int i = 0; i < estimations.size(); ++i)
	{
		if (estimations[i].size() == 0)
		{
			throw std::invalid_argument("Estimation vector should not be zero size");
		}

		if (estimations[i].size() != ground_truth[i].size())
		{
			throw std::invalid_argument("Estimation vector should be same size as a ground truth");
		}
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		residual = (estimations[i] - ground_truth[i]);
		residual = residual.array() * residual.array();
		rmse = rmse + residual;
	}

	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}

double Tools::NormalizeAngle(double &input) {
	double res = input;

	while (res >= M_PI) {
		res -= 2 * M_PI;
	}

	while (res < -M_PI) {
		res += 2 * M_PI;
	}

	return res;
}


double Tools::CalculateNIS(const VectorXd z, const VectorXd z_pred, const MatrixXd S) {
	
	double res = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
	return res;
}