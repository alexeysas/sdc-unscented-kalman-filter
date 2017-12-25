#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.5;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = M_PI_4 / 1.5;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;
	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
	is_initialized_ = false;

	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	///* Sigma point spreading parameter
	lambda_ = 3 - n_aug_;

	mes_index = 0;

	nis_radar_less_95 = 0;

	nis_radar_greater_95 = 0;

	nis_radar_95 = 7.815;

	nis_lidar_less_95 = 0;

	nis_lidar_greater_95 = 0;

	nis_lidar_95 = 5.991;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/
	/*****************************************************************************
	*  Initialization
	****************************************************************************/

	if (!is_initialized_) {
		cout << "Initialization started..." << endl;

		// if first measurement is radar - initialize from radar data by converting from polar coordinates 
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			double p = meas_package.raw_measurements_[0];
			double phi = meas_package.raw_measurements_[1];
			x_ << p * cos(phi), p * sin(phi), 0, 0, 0;

		}
		// if first measurement is laser - initialize from laser data taking  direct data
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}
		// just to ignore possible sensors we might dont know about 
		else {
			return;
		}

		// set initialization timestamp as a starting point 
		time_us_ = meas_package.timestamp_;

		// initialize object covariance matrix
		P_ << 1, 0, 0, 0, 0,
			0, 1, 0, 0, 0,
			0, 0, 1, 0, 0,
			0, 0, 0, 1, 0,
			0, 0, 0, 0, 1;

		// Initialize weights
		InitWeights();

		// done initializing, no need to predict or update
		is_initialized_ = true;

		cout << "Initialization done!" << endl;

		return;
	}

	// check if we need to switch specific sensors off
	if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) || (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_))
	{
		return;
	}

	//compute the time elapsed between the current and previous measurements
	float dt = ApplyTime(meas_package);

	// Predict
	MatrixXd Xsig_pred = Prediction(dt);

	// Update
	Update(meas_package, Xsig_pred);

	mes_index++;

	// Display NSE to console
	DisplayNIS();

}

void UKF::DisplayNIS() {
	if (mes_index % 30 == 0 && nis_radar_greater_95 > 0 && nis_radar_less_95 > 0) {
		double actual_precentage = (double)nis_radar_greater_95 / (double)(nis_radar_greater_95 + nis_radar_less_95) * 100;
		cout << "Radar NIS greater than " << nis_radar_95 << " rate:  " << actual_precentage << "%" << endl;
	}

	if (mes_index % 30 == 0 && nis_lidar_greater_95 > 0 && nis_lidar_less_95 > 0) {
		double actual_precentage = (double)nis_lidar_greater_95 / (double)(nis_lidar_greater_95 + nis_lidar_less_95) * 100;
		cout << "Lidar NIS greater than " << nis_lidar_95 << " rate:  " << actual_precentage << "%" << endl;
	}
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
MatrixXd UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	// Generate sigma points
	MatrixXd Xsig = GenerateSigmaPoints();

	// Predict sigma points
	MatrixXd  Xsig_pred = PredictSigmaPoints(Xsig, delta_t);

	// Predict mean and covariance
	PredictMeanAndCovariance(Xsig_pred);

	return Xsig_pred;
}


MatrixXd UKF::GenerateSigmaPoints() {

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.setZero();
	x_aug.head(n_x_) = x_;

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.setZero();
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_aug_ - 2, n_aug_ - 2) = std_a_ * std_a_;
	P_aug(n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;

	MatrixXd Xsig = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	MatrixXd A = P_aug.llt().matrixL();

	double sqrt_lambda = sqrt(lambda_ + n_aug_);

	//create augmented sigma points
	Xsig.col(0) = x_aug;

	for (int i = 1; i <= n_aug_; i++)
	{
		Xsig.col(i) = x_aug + A.col(i - 1) * sqrt_lambda;
		Xsig.col(i + n_aug_) = x_aug - A.col(i - 1) * sqrt_lambda;
	}

	return Xsig;
}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig, double delta_t) {

	//resulting sigma points predictions 
	MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred.setZero();

	//precalculated squared delta_t
	double delta_t_sq = delta_t * delta_t;

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd x = Xsig.col(i);
		VectorXd dx = VectorXd(n_x_);
		VectorXd noise = VectorXd(n_x_);

		double v = x(2);
		double yaw = x(3);
		double yaw_d = x(4);
		double nu_a_k = x(5);
		double nu_yawdd_k = x(6);

		noise << 0.5 * delta_t_sq * cos(yaw) * nu_a_k, 0.5 * delta_t_sq * sin(yaw) * nu_a_k,
			delta_t * nu_a_k, 0.5 * delta_t_sq * nu_yawdd_k, delta_t * nu_yawdd_k;

		if (yaw_d < 0.0001)
		{
			dx << v * cos(yaw) * delta_t, v * sin(yaw) * delta_t, 0, yaw_d * delta_t, 0;
		}
		else
		{
			dx << v / yaw_d * (sin(yaw + yaw_d * delta_t) - sin(yaw)),
				v / yaw_d * (-cos(yaw + yaw_d * delta_t) + cos(yaw)), 0, yaw_d * delta_t, 0;
		}

		Xsig_pred.col(i).head(n_x_) = x.head(n_x_) + dx + noise;
	}

	return Xsig_pred;
}


void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred) {
	x_.setZero();
	P_.setZero();
	int n_s = 2 * n_aug_ + 1;

	//predict state mean
	for (int i = 0; i < n_s; i++)
	{
		x_ += weights_(i) * Xsig_pred.col(i);
	}

	//predict state covariance matrix
	for (int i = 0; i < n_s; i++)
	{
		VectorXd residual = Xsig_pred.col(i) - x_;
		P_ += weights_(i) * residual * residual.transpose();
	}
}



void UKF::Update(MeasurementPackage meas_package, MatrixXd  Xsig_pred) {
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

		// Radar updates
		UpdateRadar(meas_package, Xsig_pred);
	}
	else {
		// Laser updates
		UpdateLidar(meas_package, Xsig_pred);
	}
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package, MatrixXd  Xsig_pred) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/

	MesPredidctions predctions = PredictLidarMeasurement(meas_package, Xsig_pred);
	UpdateStateLidar(meas_package, Xsig_pred, predctions);
}

/*
// Full version of function - replaced by optimized one below

MesPredidctions UKF::PredictLidarMeasurement(MeasurementPackage meas_package, MatrixXd  Xsig_pred) {
	//create matrix for sigma points in measurement space
	int n_z = 2;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd x = Xsig_pred.col(i);
		double px = x(0);
		double py = x(1);
		Zsig.col(i) << px, py;
	}

	//calculate mean predicted measurement
	z_pred.setZero();

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred += Zsig.col(i) * weights_(i);
	}

	//calculate covariance matrix S
	S.setZero();

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd residual = Zsig.col(i) - z_pred;
		S += residual * residual.transpose() * weights_(i);
	}

	//apply additive noise
	S(0, 0) += std_laspx_ * std_laspx_;
	S(1, 1) += std_laspy_ * std_laspy_;


	MesPredidctions result;

	result.z_pred = z_pred;
	result.S = S;
	result.Zsig = Zsig;

	cout << "Xsig" << endl;
	cout << Xsig_pred << endl;

	cout << "Zsig" << endl;
	cout << Zsig << endl;

	cout << S << endl;

	return result;
}

*/

MesPredidctions UKF::PredictLidarMeasurement(MeasurementPackage meas_package, MatrixXd  Xsig_pred) {
	//create matrix for sigma points in measurement space
	int n_z = 2;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);

	//transform sigma points into measurement space - just projection for lidar
	Zsig = Xsig_pred.block(0, 0, n_z, 2 * n_aug_ + 1);

	// calculate mean predicted measurement 
	// for lidar - no extra trasformation - so z_pred = x_pred  
	z_pred(0) = x_(0);
	z_pred(1) = x_(1);

	//calculate covariance matrix S - no need to calculate for lidar
	S = P_.block(0, 0, 2, 2);

	//apply additive noise
	S(0, 0) += std_laspx_ * std_laspx_;
	S(1, 1) += std_laspy_ * std_laspy_;

	MesPredidctions result;

	result.z_pred = z_pred;
	result.S = S;
	result.Zsig = Zsig;

	return result;
}


void UKF::UpdateStateLidar(MeasurementPackage meas_package, MatrixXd  Xsig_pred, MesPredidctions predictions) {
	int n_z = 2;

	// extract measurements
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

	// extract predictions
	VectorXd z_pred = predictions.z_pred;
	MatrixXd S = predictions.S;
	MatrixXd Zsig = predictions.Zsig;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.setZero();

	//calculate cross correlation matrix
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd res_x = Xsig_pred.col(i) - x_;
		VectorXd res_z = Zsig.col(i) - z_pred;
		Tc += res_x * res_z.transpose() * weights_(i);
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(n_x_, n_z);
	K = Tc * S.inverse();

	//update state mean and covariance matrix
	x_ = x_ + K * (z - z_pred);
	P_ = P_ - K *  S * K.transpose();

	double nis = tools.CalculateNIS(z, z_pred, S);

	if (nis > nis_lidar_95) {
		nis_lidar_greater_95++;
	}
	else {
		nis_lidar_less_95++;
	}
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package, MatrixXd  Xsig_pred) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/

	MesPredidctions predctions = PredictRadarMeasurement(meas_package, Xsig_pred);
	UpdateStateRadar(meas_package, Xsig_pred, predctions);
}

MesPredidctions UKF::PredictRadarMeasurement(MeasurementPackage meas_package, MatrixXd  Xsig_pred) {

	//create matrix for sigma points in measurement space
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd x = Xsig_pred.col(i);
		double px = x(0);
		double py = x(1);
		double v = x(2);
		double yaw = x(3);

		double ro = sqrt(px*px + py*py);
		double phi = atan2(py, px);
		double ro_dot = (px * cos(yaw) * v + py * sin(yaw)* v) / ro;
		Zsig.col(i) << ro, phi, ro_dot;
	}

	//calculate mean predicted measurement
	z_pred.setZero();

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred += Zsig.col(i) * weights_(i);
	}

	//calculate covariance matrix S
	S.setZero();

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd residual = Zsig.col(i) - z_pred;

		// normalize angle
		residual(1) = tools.NormalizeAngle(residual(1));

		S += residual * residual.transpose() * weights_(i);
	}

	//apply additive noise
	S(0, 0) += std_radr_ * std_radr_;
	S(1, 1) += std_radphi_ * std_radphi_;
	S(2, 2) += std_radrd_ * std_radrd_;

	MesPredidctions result;

	result.z_pred = z_pred;
	result.S = S;
	result.Zsig = Zsig;

	return result;
}

void UKF::UpdateStateRadar(MeasurementPackage meas_package, MatrixXd  Xsig_pred, MesPredidctions predictions) {

	int n_z = 3;

	// extract measurements
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

	// extract predictions
	VectorXd z_pred = predictions.z_pred;
	MatrixXd S = predictions.S;
	MatrixXd Zsig = predictions.Zsig;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.setZero();

	//calculate cross correlation matrix
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd res_x = Xsig_pred.col(i) - x_;
		VectorXd res_z = Zsig.col(i) - z_pred;

		res_x(3) = tools.NormalizeAngle(res_x(3));
		res_z(1) = tools.NormalizeAngle(res_z(1));

		Tc += res_x * res_z.transpose() * weights_(i);
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(n_x_, n_z);
	K = Tc * S.inverse();

	//update state mean and covariance matrix
	x_ = x_ + K * (z - z_pred);
	P_ = P_ - K *  S * K.transpose();

	//calculatre nis and rate when it is greater than 95% threashold
	double nis = tools.CalculateNIS(z, z_pred, S);

	if (nis > nis_radar_95) {
		nis_radar_greater_95++;
	}
	else {
		nis_radar_less_95++;
	}
}




double UKF::ApplyTime(MeasurementPackage meas_package) {
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

	return dt;
}

void UKF::InitWeights() {

	cout << "Initializing weights...." << endl;
	int n_s = 2 * n_aug_ + 1;
	weights_ = VectorXd(n_s);

	//set weights
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i < n_s; i++)
	{
		weights_(i) = 1 / (2 * (lambda_ + n_aug_));
	}
}