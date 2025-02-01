//
// Created by VS0121 on 19/8/2024.
//

#include "Cylinder.h"

double Cylinder::F(const double x_, const double y_) const
{
    return x_*x_ + y_*y_ - 0.25 * pow(D,2);
}

Eigen::Vector2d Cylinder::GradF(const double x_, const double y_) const
{
    Eigen::Vector2d grad = {2*x_, 2*y_};
    return grad;
}

double Cylinder::DistanceFromWall(const double x0_, const double y0_) const {
    Eigen::Vector2d x0(x0_,y0_);
    return x0.norm() - D/2;
}
