//
// Created by VS0121 on 15/9/2024.
//

#include "CalculateFields.h"

int sum(int n)
{
    int SUM = 0;
    for (int i = 0; i <= n; ++i) {
        SUM += i;
    }
    return SUM;
}

// Function to load (x, y, T) data from all particles at a given timestep
Eigen::MatrixXd loadTemperatureDataAtTimestep(const std::vector<std::string>& filenames, int timestep) {
    std::vector<Eigen::Vector3d> data;  // Store x, y, T for all particles at the timestep

    for (const auto& filename : filenames) {
        std::ifstream file(filename);
        std::string line;

        // Skip header
        std::getline(file, line);

        // Read each line until we reach the desired timestep
        for (int i = 0; i <= timestep; ++i) {
            std::getline(file, line);  // Read timestep line
        }

        //std::cout << line;

        // Parse the line for the current timestep
        std::stringstream ss(line);
        double t, x, y, z, u, v, m, T;
        ss >> t >> x >> y >> z >> u >> v >> m >> T;

        // Store (x, y, T)
        data.push_back(Eigen::Vector3d(x, y, T));
    }

    // Convert to Eigen matrix
    Eigen::MatrixXd result(data.size(), 3);
    for (size_t i = 0; i < data.size(); ++i) {
        result.row(i) = data[i];
    }

    return result;
}

// Function to fit a polynomial to the data using least squares
Eigen::VectorXd fitPolynomial(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values, int order) {
    int M = points.rows();  // Number of data points
    int num_coeffs = (order + 1) * (order + 2) / 2;  // Number of coefficients for 2D polynomial

    Eigen::MatrixXd A(M, num_coeffs);  // Design matrix

    // Build design matrix for polynomial fitting
//    int col = 0;
//    for (int i = 0; i <= order; ++i) {
//        for (int j = 0; j <= order - i; ++j) {
//            for (int k = 0; k < M; ++k) {
//                A(k, col) = std::pow(points(k, 0), i) * std::pow(points(k, 1), j);
//            }
//            ++col;
//        }
//    }
    auto x_col = points.col(0);
    auto y_col = points.col(1);

    int j = 0;
    for (int n = 0; n <= order; ++n) {
        for (int jn = 0; jn <= n; ++jn) {
                //int j = sum(n) + jn;
                int exp_x = n-jn;
                int exp_y = jn;
                A.col(j) = x_col.array().pow(exp_x) * y_col.array().pow(exp_y);
                j++;
        }
    }



    // Solve for polynomial coefficients using least squares
    Eigen::VectorXd coeffs = A.colPivHouseholderQr().solve(T_values);
    return coeffs;
}



// Function to evaluate the fitted polynomial at a given (x, y)
double evaluatePolynomial(const Eigen::VectorXd& coeffs, double x, double y, int order) {
    double value = 0.0;
//    int col = 0;
//    for (int i = 0; i <= order; ++i) {
//        for (int j = 0; j <= order - i; ++j) {
//            value += coeffs(col) * std::pow(x, i) * std::pow(y, j);
//            ++col;
//        }
//    }
    int j = 0;
    for (int n = 0; n <= order; ++n) {
        for (int jn = 0; jn <= n; ++jn) {
            value += pow(x, n-jn) * pow(y, jn);
            j++;
        }
    }
    return value;
}

// Function to calculate RMSE between fitted polynomial and original data
double calculateRMSE(const Eigen::VectorXd& T_values, const Eigen::VectorXd& predicted_T) {
    Eigen::VectorXd diff = T_values - predicted_T;
    //std::cout << "T = " << T_values.transpose() << "\n";
    //std::cout << "T_ = " << predicted_T.transpose() << "\n";

    return std::sqrt(diff.squaredNorm() / T_values.size());
}

// Function to find the best polynomial order (between 2 and 10)
int findBestPolynomialOrder(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values, int min_order, int max_order, Eigen::VectorXd& best_coeffs) {
    int best_order = min_order;
    double best_rmse = std::numeric_limits<double>::max();

    // Loop over orders from min_order to max_order
    for (int order = min_order; order <= max_order; ++order) {
        Eigen::VectorXd coeffs = fitPolynomial(points, T_values, order);

        // Predict temperatures using the fitted polynomial
        Eigen::VectorXd predicted_T(T_values.size());
        for (int i = 0; i < points.rows(); ++i) {
            predicted_T(i) = evaluatePolynomial(coeffs, points(i, 0), points(i, 1), order);
        }

        // Calculate RMSE for the current order
        double rmse = calculateRMSE(T_values, predicted_T);

        // Update best order if the current one is better
        if (rmse < best_rmse) {
            best_rmse = rmse;
            best_order = order;
            best_coeffs = coeffs;
        }
    }

    return best_order;
}

// Function to calculate the mean temperature field from polynomial fits
Eigen::MatrixXd calculateMeanTemperatureField(const std::vector<Eigen::VectorXd>& coeffs_list,
                                              const std::vector<int>& best_orders,
                                              int grid_size,
                                              const double x_min,
                                              const double x_max,
                                              const double y_min,
                                              const double y_max)
{
    Eigen::MatrixXd mean_field(grid_size * grid_size, 3);
    int idx = 0;
    double dx = (x_max - x_min) / (grid_size - 1);
    double dy = (y_max - y_min) / (grid_size - 1);

    // Loop over grid points and compute the mean temperature
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            double T_sum = 0.0;

            // Sum temperature values from all polynomial fits
            for (size_t t = 0; t < coeffs_list.size(); ++t) {
                T_sum += evaluatePolynomial(coeffs_list[t], x, y, best_orders[t]);
            }

            double T_mean = T_sum / coeffs_list.size();
            mean_field(idx, 0) = x;
            mean_field(idx, 1) = y;
            mean_field(idx, 2) = T_mean;
            ++idx;
        }
    }
    return mean_field;
}

// Function to save the mean temperature field for Gnuplot
void saveTemperatureField(const Eigen::MatrixXd& field, const std::string& filename) {
    std::ofstream file(filename);
    for (int i = 0; i < field.rows(); ++i) {
        file << field(i, 0) << " " << field(i, 1) << " " << field(i, 2) << "\n";
    }
}

// 2D Gaussian model function
double gaussian2D(double A, double x0, double y0, double sigma_x, double sigma_y, double x, double y) {
    return A * std::exp(-((x - x0) * (x - x0) / (2 * sigma_x * sigma_x) + (y - y0) * (y - y0) / (2 * sigma_y * sigma_y)));
}

// Function to fit a 2D Gaussian to the data using nonlinear least squares
Eigen::VectorXd fitGaussian2D(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values) {
    int M = points.rows();  // Number of data points

    // Initial guesses for the parameters [A, x0, y0, sigma_x, sigma_y]
    double A = T_values.maxCoeff();  // Amplitude starts as the max temperature
    double x0 = points.col(0).mean();  // Initial guess for the center x
    double y0 = points.col(1).mean();  // Initial guess for the center y
    double sigma_x = 1.0;
    double sigma_y = 1.0;

    // Perform nonlinear least squares (simple iterative method)
    double learning_rate = 0.01;  // Learning rate for gradient descent
    int max_iters = 1000;  // Maximum number of iterations
    for (int iter = 0; iter < max_iters; ++iter) {
        Eigen::VectorXd residuals(M);
        Eigen::MatrixXd J(M, 5);  // Jacobian matrix (M data points, 5 parameters)

        for (int i = 0; i < M; ++i) {
            double x = points(i, 0);
            double y = points(i, 1);
            double G = gaussian2D(A, x0, y0, sigma_x, sigma_y, x, y);

            residuals(i) = T_values(i) - G;

            // Calculate partial derivatives for the Jacobian
            J(i, 0) = -G / A;  // dG/dA
            J(i, 1) = G * (x - x0) / (sigma_x * sigma_x);  // dG/dx0
            J(i, 2) = G * (y - y0) / (sigma_y * sigma_y);  // dG/dy0
            J(i, 3) = G * (x - x0) * (x - x0) / (sigma_x * sigma_x * sigma_x);  // dG/dsigma_x
            J(i, 4) = G * (y - y0) * (y - y0) / (sigma_y * sigma_y * sigma_y);  // dG/dsigma_y
        }

        // Update the parameters using gradient descent
        Eigen::VectorXd delta = (J.transpose() * J).ldlt().solve(J.transpose() * residuals);
        A -= learning_rate * delta(0);
        x0 -= learning_rate * delta(1);
        y0 -= learning_rate * delta(2);
        sigma_x -= learning_rate * delta(3);
        sigma_y -= learning_rate * delta(4);

        // Stopping condition
        if (delta.norm() < 1e-6) {
            break;
        }
    }

    Eigen::VectorXd params(5);
    params << A, x0, y0, sigma_x, sigma_y;
    return params;
}

// Function to evaluate the fitted 2D Gaussian at a given (x, y)
double evaluateGaussian2D(const Eigen::VectorXd& params, double x, double y) {
    double A = params(0);
    double x0 = params(1);
    double y0 = params(2);
    double sigma_x = params(3);
    double sigma_y = params(4);
    return gaussian2D(A, x0, y0, sigma_x, sigma_y, x, y);
}

// Function to find the best 2D Gaussian fit (only one Gaussian)
Eigen::VectorXd findBestGaussianFit(const Eigen::MatrixXd& points, const Eigen::VectorXd& T_values) {
    Eigen::VectorXd best_params = fitGaussian2D(points, T_values);

    // Predict temperatures using the fitted Gaussian
    Eigen::VectorXd predicted_T(T_values.size());
    for (int i = 0; i < points.rows(); ++i) {
        predicted_T(i) = evaluateGaussian2D(best_params, points(i, 0), points(i, 1));
    }

    // Calculate RMSE for the Gaussian fit
    double rmse = calculateRMSE(T_values, predicted_T);

    std::cout << "Gaussian Fit RMSE: " << rmse << std::endl;
    return best_params;
}

// Function to calculate the mean temperature field from Gaussian fits
Eigen::MatrixXd calculateMeanTemperatureFieldGaussian(const std::vector<Eigen::VectorXd>& params_list,
                                              int grid_size,
                                              const double x_min,
                                              const double x_max,
                                              const double y_min,
                                              const double y_max) {
    Eigen::MatrixXd mean_field(grid_size * grid_size, 3);
    int idx = 0;
    double dx = (x_max - x_min) / (grid_size - 1);
    double dy = (y_max - y_min) / (grid_size - 1);

    // Loop over grid points and compute the mean temperature
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            double T_sum = 0.0;

            // Sum temperature values from all Gaussian fits
            for (size_t t = 0; t < params_list.size(); ++t) {
                T_sum += evaluateGaussian2D(params_list[t], x, y);
            }

            double T_mean = T_sum / params_list.size();
            mean_field(idx, 0) = x;
            mean_field(idx, 1) = y;
            mean_field(idx, 2) = T_mean;
            ++idx;
        }
    }
    return mean_field;
}
