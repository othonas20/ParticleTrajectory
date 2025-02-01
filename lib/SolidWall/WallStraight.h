//
// Created by VS0121 on 23/8/2024.
//

#ifndef PARTICLETRAJECTORY_WALLSTRAIGHT_H
#define PARTICLETRAJECTORY_WALLSTRAIGHT_H

#include "SolidWall.h"

class WallStraight : public SolidWall{
public:
    WallStraight(const Eigen::Vector2d pointSolid_, const double delta_critical_)
    :SolidWall(pointSolid_,delta_critical_){}

    virtual double F(const double x_, const double y_) const;

    virtual Eigen::Vector2d GradF(const double x_, const double y_) const;

    virtual double DistanceFromWall(const double x0_, const double y0_) const;

    virtual ~WallStraight() = default;

};


#endif //PARTICLETRAJECTORY_WALLSTRAIGHT_H
