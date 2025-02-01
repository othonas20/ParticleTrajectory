//
// Created by VS0121 on 16/1/2024.
//

#ifndef SOLIDWALL_SOLIDWALL_H
#define SOLIDWALL_SOLIDWALL_H

#include <iostream>
#include <Eigen/Dense>

// abstract class for solid walls
// double F and Vector2d gradF functions need to be defined
class SolidWall
{
private:
    // 2D point inside solid region
    const Eigen::Vector2d pointSolid;
    const double delta_critical;

public:

    SolidWall(const Eigen::Vector2d pointSolid_, const double delta_critical_);

    virtual double F(const double x, const double y) const = 0;

    virtual Eigen::Vector2d GradF(const double x, const double y) const = 0;

    virtual double DistanceFromWall(const double x0_, const double y0_) const = 0;

    Eigen::Vector2d Normal(const double x, const double y) const;


    Eigen::Vector2d Tangent(const double x, const double y) const;


    double Fcrit(const double x, const double y) const;


    Eigen::Vector2d getPointSolid() const;

    double getDeltaCritical() const;

    bool PointInFluidRegion(const double x, const double y) const { return F(x,y) > 0;}

    virtual ~SolidWall() = default;
};

#endif //SOLIDWALL_SOLIDWALL_H
