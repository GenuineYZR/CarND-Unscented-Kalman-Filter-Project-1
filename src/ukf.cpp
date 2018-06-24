#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
	// initially set to false, set to true in first call of ProcessMeasurement
	is_initialized_ = false;
  
	// if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 1, 0,
	  0, 0, 0, 0, 1;

  // initial timestamp
  time_us_ = 0.0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);

  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (int i = 1; i < 2 * n_aug_ + 1; i++)
  {
	  weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
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
	if (!is_initialized_)
	{
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			double a = meas_package.raw_measurements_[0];
			double b = meas_package.raw_measurements_[1];
			double py = a * sin(b);
			double px = a * cos(b);
			x_ << px, py, 0, 0, 0;
		}
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}

	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	
	Prediction(delta_t);
	
	if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	{
		UpdateLidar(meas_package);
	}

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		UpdateRadar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	VectorXd x_aug = VectorXd(7);
	x_aug.fill(0.0);

	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	MatrixXd P_aug = MatrixXd(7, 7);
	P_aug.fill(0.0);

	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	MatrixXd x_sig = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	x_sig.fill(0.0);

	MatrixXd A = P_aug.llt().matrixL();

	x_sig.col(0) = x_aug;

	for (int i = 0; i < n_aug_; i++)
	{
		x_sig.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		x_sig.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}
	
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred_.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd x = VectorXd(7);
		x.fill(0.0);

		x = x_sig.col(i);

		double px = x(0);
		double py = x(1);
		double v = x(2);
		double yaw = x(3);
		double yawd = x(4);
		double nu_a_ = x(5);
		double nu_yawdd_ = x(6);
		
		if (fabs(yawd) > 0.0001)
		{
			px = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw)) + 0.5 * delta_t * delta_t * cos(yaw) * nu_a_;
			py = py + v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw)) + 0.5 * delta_t * delta_t * sin(yaw) * nu_a_;
			v = v + delta_t * nu_a_;
			yaw = yaw + yawd * delta_t + 0.5 * nu_yawdd_ * delta_t * delta_t;
			yawd = yawd + delta_t * nu_yawdd_;
		}
		else
		{
			px = px + v * cos(yaw) * delta_t + 0.5 * delta_t * delta_t * cos(yaw) * nu_a_;
			py = py + v * sin(yaw) * delta_t + 0.5 * delta_t * delta_t * sin(yaw) * nu_a_;
			v = v + delta_t * nu_a_;
			yaw = yaw + yawd * delta_t + 0.5 * nu_yawdd_ * delta_t * delta_t;
			yawd = yawd + delta_t * nu_yawdd_;
		}

		Xsig_pred_.col(i) << px, py, v, yaw, yawd;
	}

	x_.fill(0.0);
	P_.fill(0.0);
	//predicted state mean
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	//predicted state covariance matrix
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		// state difference
		VectorXd x_diff = VectorXd(2 * n_aug_ + 1);
		x_diff.fill(0.0);

		x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		if (x_diff(3) > M_PI)
		{
			x_diff(3) = x_diff(3) - 2 * M_PI;
		}
		if (x_diff(3) < - M_PI)
		{
			x_diff(3) = x_diff(3) + 2 * M_PI;
		}
		
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //set measurement dimension, lidar can measure px and py
	int n_z = 2;
	
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	
	// set measurement matrix
	MatrixXd H = MatrixXd(n_z, n_x_);
	H << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0;
	
	z_pred = H * x_;
	
	// set measurement covariance matrix
	S = H * P_ * H.transpose();
	
	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_ * std_laspy_;
	
	S = S + R;
	
	//Kalman gain K
	MatrixXd K = P_ * H.transpose() * S.inverse();
	
	// update state vector
	VectorXd z = VectorXd(n_z);
	double a = meas_package.raw_measurements_[0];
	double b = meas_package.raw_measurements_[1];
	
	z << a, b;
	x_ = x_ + K * (z - z_pred);
	
	// upate covariance matrix
	MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
	
	P_ = (I - K * H) * P_;
	
	// calculate the lidar NIS
	double NIS = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;
	
  //create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0.0);
	
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd Z = VectorXd(n_x_);

		Z = Xsig_pred_.col(i);

		double px = Z(0);
		double py = Z(1);
		double v = Z(2);
		double yaw = Z(3);
		double yawd = Z(4);

		double r = sqrt(px * px + py * py);
		double phi, rd;

		if (fabs(px) < 0.0001)
		{
			cout << "When converting x_ from polar to cartesian coordinates - Division by Zero" << endl;
		}
		else
		{
			phi = atan2(py, px);
		}

		if (fabs(r) < 0.0001)
		{
			cout << "When converting x_ from polar to cartesian coordinates - Division by Zero" << endl;
		}
		else
		{
			rd = (px * v * cos(yaw) + py * v * sin(yaw)) / r;
		}

		Zsig.col(i) << r, phi, rd;
	}
	

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//innovation covariance matrix S
	S.fill(0.0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		
		//angle normalization
		while (z_diff(1)> M_PI)
		{
			z_diff(1) -= 2.*M_PI;
		}
		while (z_diff(1)<-M_PI)
		{
			z_diff(1) += 2.*M_PI;
		}

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}
	
	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R.fill(0.0);
	R(0, 0) = std_radr_ * std_radr_;
	R(1, 1) = std_radphi_ * std_radphi_;
	R(2, 2) = std_radrd_ * std_radrd_;

	S = S + R;
	
	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		// state difference
		VectorXd x_diff = VectorXd(2 * n_aug_ + 1);
		x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		if (x_diff(3) > M_PI)
		{
			x_diff(3) = x_diff(3) - 2 * M_PI;
		}
		if (x_diff(3) < -M_PI)
		{
			x_diff(3) = x_diff(3) + 2 * M_PI;
		}
		
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI)
		{
			z_diff(1) -= 2.*M_PI;
		}
		while (z_diff(1)<-M_PI)
		{
			z_diff(1) += 2.*M_PI;
		}

		Tc = Tc + x_diff * z_diff.transpose() * weights_(i);
	}
	
	//Kalman gain K
	MatrixXd K = Tc * S.inverse();
	
	// update state vector
	VectorXd z = VectorXd(n_z);
	
	double a = meas_package.raw_measurements_[0];
	
	double b = meas_package.raw_measurements_[1];
	
	double c = meas_package.raw_measurements_[2];
	
	z << a, b, c;
	x_ = x_ + K * (z - z_pred);
	
	// upate covariance matrix
	P_ = P_ - K * S * K.transpose();
	
	// calculate the lidar NIS
	double NIS = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}
