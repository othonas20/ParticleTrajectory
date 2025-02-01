//
// Created by VS0121 on 23/8/2024.
//

#include "WallStraight.h"

double WallStraight::F(const double x_, const double y_) const {
    return y_;
}

Eigen::Vector2d WallStraight::GradF(const double x_, const double y_) const {
    return Eigen::Vector2d(0,1);
}

double WallStraight::DistanceFromWall(const double x0_, const double y0_) const {
    return y0_;
}
