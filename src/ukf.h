#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct MesPredidctions {
public:
	MatrixXd  Zsig;
	VectorXd  z_pred; 
	MatrixXd  S;
};

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  // Counter for the measurements
  int mes_index;

  // Counter for Radar NIS greater than 7.8
  int nis_radar_less_95;

  // Counter for Radar NIS less than 7.8
  int nis_radar_greater_95;

  // 7.8
  double nis_radar_95;

  // Counter for Lidar NIS less than 5.9
  int nis_lidar_less_95;

  // Counter for Lidar NIS greater than 5.9
  int nis_lidar_greater_95;

  // 5.9
  double nis_lidar_95;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  MatrixXd Prediction(double delta_t);

  MatrixXd GenerateSigmaPoints();

  MatrixXd PredictSigmaPoints(MatrixXd Xsig, double delta_t);

  void PredictMeanAndCovariance(MatrixXd Xsig_pred);

  void Update(MeasurementPackage meas_package, MatrixXd  Xsig_pred);

  void DisplayNIS();

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package, MatrixXd  Xsig_pred);

  MesPredidctions PredictLidarMeasurement(MeasurementPackage meas_package, MatrixXd  Xsig_pred);

  void UpdateStateLidar(MeasurementPackage meas_package, MatrixXd  Xsig_pred, MesPredidctions predictions);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package, MatrixXd  Xsig_pred);

  MesPredidctions PredictRadarMeasurement(MeasurementPackage meas_package, MatrixXd  Xsig_pred);

  void UpdateStateRadar(MeasurementPackage meas_package, MatrixXd  Xsig_pred, MesPredidctions predictions);

  double ApplyTime(MeasurementPackage meas_package);

  void InitWeights();

private:
	Tools tools;
};

#endif /* UKF_H */
