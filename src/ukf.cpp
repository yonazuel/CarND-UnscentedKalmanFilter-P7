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
 */
UKF::UKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.75;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Number of sigma points
  n_sig_ = 15;

  // intial sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // sigma points parameter
  lambda_ = 3 - n_aug_;

  // weights of sigma points
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  double w_i = 0.5 / (lambda_ + n_aug_);
  for(int i = 1; i < n_sig_; i++) {
    weights_(i) = w_i;
  }

  // initial radar measurement covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  // initial lidar measurement covariance matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  time_us_ = 0;
}

UKF::~UKF() {}

double UKF::Mod2PI(double ang) {
    double ang_norm = ang;
    while (ang_norm > M_PI) ang_norm -= 2. * M_PI;
    while (ang_norm < -M_PI) ang_norm += 2. * M_PI;
    return ang_norm;
}

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

  /* if it is the first measurement, I use it to initialize the state
  vector and covariance matrix */
  if(!is_initialized_) {

    // initial covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 0.5, 0, 0,
          0, 0, 0, 0.5, 0,
          0, 0, 0, 0, 1;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      //Convert radar from polar to cartesian coordinates and initialize state.
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rhod = meas_package.raw_measurements_[2];
      
      double px = rho * cos(phi);
      double py = rho * sin(phi);

      x_ << px, py, rhod, 0, 0;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    if(fabs(x_[0]) < 0.0001) {
      x_[0] = 0.0001;
    }
    if(fabs(x_[1]) < 0.0001) {
      x_[1] = 0.0001;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Here I compute the timestep and go through the Prediction and Update parts
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
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

  // First I create the augmented state, covariance matrix
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  // Then I compute the augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;
  MatrixXd L = P_aug.llt().matrixL();
  double slna = sqrt(lambda_ + n_aug_);
  VectorXd Col_sig_aug;
  for(int i = 0; i < n_aug_; i++){
    Col_sig_aug = slna * L.col(i);
    Xsig_aug.col(i+1) = x_aug + Col_sig_aug;
    Xsig_aug.col(i+1+n_aug_) = x_aug - Col_sig_aug;
  }

  // Now I need to predict the sigma points
  for(int i = 0; i < n_sig_; i++) {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double px_p, py_p, v_p, yaw_p, yawd_p;

    double cos_yaw = cos(yaw);
    double sin_yaw = sin(yaw);
    double dt2 = delta_t * delta_t;

    yaw_p = yaw + yawd * delta_t;

    if(fabs(yawd) > 0.0001) {
      double vyawd = v/yawd;
      px_p = p_x + vyawd * (sin(yaw_p) - sin_yaw);
      py_p = p_y + vyawd * (-cos(yaw_p) + cos_yaw);
    }
    else {
      double vdt = v * delta_t;
      px_p = p_x + vdt * cos_yaw;
      py_p = p_y + vdt * sin_yaw;
    }

    v_p = v;
    yawd_p = yawd;

    px_p += 0.5 * dt2 * nu_a * cos_yaw;
    py_p += 0.5 * dt2 * nu_a * sin_yaw;
    v_p += delta_t * nu_a;
    yaw_p += 0.5 * dt2 * nu_yawdd;
    yawd_p += delta_t * nu_yawdd;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // So the predicted state mean is
  x_ = Xsig_pred_ * weights_;
  // The predicted state covariance is
  P_.fill(0);
  for(int i = 0; i < n_sig_; i++) {
    VectorXd x_dif = Xsig_pred_.col(i) - x_;
    x_dif(3) = Mod2PI(x_dif(3));
    P_ = P_ + weights_(i) * x_dif * x_dif.transpose();
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

  // I create the measurement sigma points
  MatrixXd Zsig = Xsig_pred_.block(0, 0, 2, n_sig_);
  UpdateUKF(meas_package, Zsig, 2);
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

  // I create the measurement sigma points
  MatrixXd Zsig = MatrixXd(3, n_sig_);
  for(int i = 0; i < n_sig_; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    
    Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) = (p_x*v*cos(yaw) + p_y*v*sin(yaw)) / Zsig(0,i);
  }

  UpdateUKF(meas_package, Zsig, 3);
}

/** 
 * Updates the state and covariance matrix
 */
void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z) {
  // First I predict the measurement mean
  VectorXd z_pred = VectorXd(n_z);
  z_pred = Zsig * weights_;

  // Then I predict the measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0);
  for(int i = 0; i < n_sig_; i++) {
    VectorXd z_dif = Zsig.col(i) - z_pred;
    z_dif(1) = Mod2PI(z_dif(1));
    S = S + weights_(i) * z_dif * z_dif.transpose();
  }

  // The measurement covariance noise is
  MatrixXd R = MatrixXd(n_z, n_z);
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    R = R_radar_;
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    R = R_lidar_;
  }
  S = S + R;

  // Now I compute the cross correlation matrix
  MatrixXd T = MatrixXd(n_x_, n_z);
  T.fill(0);
  for(int i = 0; i < n_sig_; i++) {
    VectorXd x_dif = Xsig_pred_.col(i) - x_;
    x_dif(3) = Mod2PI(x_dif(3));
    VectorXd z_dif = Zsig.col(i) - z_pred;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      z_dif(1) = Mod2PI(z_dif(1));
    }
    T = T + weights_(i) * x_dif * z_dif.transpose();
  }

  // The Kalman gain is
  MatrixXd K = T * S.inverse();

  // Finally I update the state and covariance matrix
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_dif = z - z_pred;
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    z_dif(1) = Mod2PI(z_dif(1));
  }

  x_ = x_ + K * z_dif;
  P_ = P_ - K * S * K.transpose();

  // I also compute the NIS
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    nis_radar_ = z_dif.transpose() * S.inverse() * z_dif;
    //cout<<nis_radar_<<endl;
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    nis_laser_ = z_dif.transpose() * S.inverse() * z_dif;
    //cout<<nis_laser_<<endl;
  }
}
