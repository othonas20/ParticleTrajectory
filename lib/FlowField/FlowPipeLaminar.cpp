//
// Created by VS0121 on 23/8/2024.
//

#include "FlowPipeLaminar.h"

double FlowPipeLaminar::u(const double x_, const double y_) const {
    return (PointInsideBounds(x_,y_)) ? Umax * (1 - pow(2*(y_-D/2)/D, 2)) : NAN;
}

double FlowPipeLaminar::v(const double x_, const double y_) const {
    return (PointInsideBounds(x_,y_)) ? 0 : NAN;
}

Eigen::Matrix2d FlowPipeLaminar::dVjdxi(const double x_, const double y_) const {
    if (!PointInsideBounds(x_,y_))
        return Eigen::Matrix2d::Ones() * NAN;
    else {
        Eigen::Matrix2d gradV = Eigen::Matrix2d::Zero();
        gradV(1, 0) = -4 * Umax * (2*y_-D) / pow(D, 2); // du/dy
        return gradV;
    }
}

double FlowPipeLaminar::p(const double x_, const double y_) const
{
    const double l = x_ - this->xMin;
    const double u_mean = this->Umax / 3;
    return this->p0 - 32 * this->mu * l * u_mean / pow(this->D, 2);
}

bool FlowPipeLaminar::PointInsideBounds(const double x_, const double y_) const {
    return y_ >= 0 && y_ <= D;
}
