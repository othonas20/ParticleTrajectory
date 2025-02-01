//
// Created by VS0121 on 23/8/2024.
//

#include "Pipe.h"

double Pipe::F(const double x_, const double y_) const {
    return y_*y_; //- 0.25*D*D;
}

Eigen::Vector2d Pipe::GradF(const double x_, const double y_) const {
    return Eigen::Vector2d(0, 2*y_);
}

double Pipe::DistanceFromWall(const double x0_, const double y0_) const {
    return abs(std::min(y0_, D - y0_));
}
