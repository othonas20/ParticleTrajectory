//
// Created by VS0121 on 15/9/2024.
//

#ifndef PARTICLETRAJECTORY_CALCULATEFIELDS_H
#define PARTICLETRAJECTORY_CALCULATEFIELDS_H

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <string>
#include <cmath>

// Function to load (x, y, T) data from all particles at a given timestep
Eigen::MatrixXd loadTemperatureDataAtTimestep(const std::vector<std::string>& filenames, int timestep) ;

// Function to fit a polynomial to the data using least squares
Eigen::VectorXd fitPolynomial(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values, int order);

// Function to evaluate the fitted polynomial at a given (x, y)
double evaluatePolynomial(const Eigen::VectorXd& coeffs, double x, double y, int order);

// Function to calculate RMSE between fitted polynomial and original data
double calculateRMSE(const Eigen::VectorXd& T_values, const Eigen::VectorXd& predicted_T);

// Function to find the best polynomial order (between 2 and 10)
int findBestPolynomialOrder(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values, int min_order, int max_order, Eigen::VectorXd& best_coeffs);

// Function to calculate the mean temperature field from polynomial fits
Eigen::MatrixXd calculateMeanTemperatureField(const std::vector<Eigen::VectorXd>& coeffs_list,
                                              const std::vector<int>& best_orders,
                                              int grid_size,
                                              const double xMin,
                                              const double xMax,
                                              const double yMin,
                                              const double yMax);

// Function to save the mean temperature field for Gnuplot
void saveTemperatureField(const Eigen::MatrixXd& field, const std::string& filename);

// 2D Gaussian model function
double gaussian2D(double A, double x0, double y0, double sigma_x, double sigma_y, double x, double y);
// Function to fit a 2D Gaussian to the data using nonlinear least squares
Eigen::VectorXd fitGaussian2D(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values);

// Function to evaluate the fitted 2D Gaussian at a given (x, y)
double evaluateGaussian2D(const Eigen::VectorXd& params, double x, double y);

// Function to find the best 2D Gaussian fit (only one Gaussian)
Eigen::VectorXd findBestGaussianFit(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values);

// Function to calculate the mean temperature field from Gaussian fits
Eigen::MatrixXd calculateMeanTemperatureFieldGaussian(const std::vector<Eigen::VectorXd>& params_list,
                                              int grid_size,
                                              const double x_min,
                                              const double x_max,
                                              const double y_min,
                                              const double y_max);
 

#endif //PARTICLETRAJECTORY_CALCULATEFIELDS_H
