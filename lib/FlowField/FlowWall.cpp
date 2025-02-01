//
// Created by VS0121 on 20/8/2024.
//

#include "FlowWall.h"

double FlowWall::Psi(const double x_, const double y_) const
{
    return (y_ > 0) ? 2 * A * x_ * y_ : 0;
}

double FlowWall::u(const double x_, const double y_) const
{
    return (y_ > 0) ? 2 * A * x_ : 0;
}

double FlowWall::v(const double x_, const double y_) const
{
    return (y_ > 0) ? - 2 * A * y_ : 0;
}

Eigen::Matrix2d FlowWall::dVjdxi(const double x_, const double y_) const
{
    if (y_ <= 0)
        return Eigen::Matrix2d::Zero();

    const double e = 1e-4;
    Eigen::Matrix2d gradV;
    gradV(0,0) = (u(x_ + e, y_) - u(x_ - e, y_)) / (2*e); // du/dx
    gradV(0,1) = (v(x_ + e, y_) - v(x_ - e, y_)) / (2*e); // dv//dx
    gradV(1,0) = (u(x_, y_ + e) - u(x_, y_ - e))/ (2*e); // du/dy
    gradV(1,1) = (v(x_, y_ + e) - v(x_, y_ - e)) / (2*e); // dv/dy

    return gradV;
}