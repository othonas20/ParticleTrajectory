//
// Created by VS0121 on 20/8/2024.
//

#ifndef PARTICLETRAJECTORY_FLOWWALL_H
#define PARTICLETRAJECTORY_FLOWWALL_H

#include "FlowField.h"

class FlowWall : public FlowField{
private:
    double A;

public:

    FlowWall(const double A_,
             const double xMin_,
             const double xMax_,
             const double yMin_,
             const double yMax_,
             const double rho_,
             const double p0_,
             const double mu_,
             const double U0_
    ):FlowField(xMin_,xMax_,yMin_,yMax_,rho_,p0_,mu_,U0_),
      A(A_){}

    virtual double Psi(const double x_, const double y_) const;

    virtual double u(const double x_, const double y_) const;

    virtual double v(const double x_, const double y_) const;

    virtual Eigen::Matrix2d dVjdxi(const double x_, const double y_) const;

    virtual ~FlowWall() = default;

};


#endif //PARTICLETRAJECTORY_FLOWWALL_H
