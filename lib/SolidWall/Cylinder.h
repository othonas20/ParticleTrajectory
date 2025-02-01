//
// Created by VS0121 on 19/8/2024.
//

#ifndef PARTICLETRAJECTORY_CYLINDER_H
#define PARTICLETRAJECTORY_CYLINDER_H

#include "SolidWall.h"

class Cylinder : public SolidWall{
private:
    double D;

public:
    Cylinder(const Eigen::Vector2d pointSolid_, const double delta_critical_, const double D_)
            : SolidWall(pointSolid_,delta_critical_), D(D_){}

    virtual double F(const double x_, const double y_) const;

    virtual Eigen::Vector2d GradF(const double x_, const double y_) const;

    virtual double DistanceFromWall(const double x0_, const double y0_) const;

    virtual ~Cylinder() = default;
};


#endif //PARTICLETRAJECTORY_CYLINDER_H
