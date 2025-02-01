//
// Created by VS0121 on 19/1/2024.
//

#include "SolidWall.h"

SolidWall::SolidWall(const Eigen::Vector2d pointSolid_, const double delta_critical_):delta_critical(delta_critical_)
{
    std::cout << "SolidWall object constructed!\n"
                 "Dont forget to fill solid line F(x,y)!\n";
}

Eigen::Vector2d SolidWall::Normal(const double x, const double y) const
{
    Eigen::Vector2d gradF_ = GradF(x,y);
    return gradF_ / gradF_.norm();
}

Eigen::Vector2d SolidWall::Tangent(const double x, const double y) const
{
    Eigen::Vector2d gradF_ = GradF(x,y);

    Eigen::Vector2d returnVec(1, -gradF_.x() / gradF_.y());
    returnVec /= returnVec.norm();
    return returnVec;

}

double SolidWall::Fcrit(const double x, const double y) const
{
    Eigen::Vector2d gradF_ = GradF(x,y);
    return delta_critical * double(gradF_.norm());
}

Eigen::Vector2d SolidWall::getPointSolid() const { return pointSolid;}

double SolidWall::getDeltaCritical() const { return delta_critical;}