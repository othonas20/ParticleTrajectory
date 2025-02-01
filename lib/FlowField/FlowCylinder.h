//
// Created by VS0121 on 19/8/2024.
//

#ifndef PARTICLETRAJECTORY_FLOWCYLINDER_H
#define PARTICLETRAJECTORY_FLOWCYLINDER_H

#include "FlowField.h"

class FlowCylinder : public FlowField{
private:
    double D;
    double xc;
    double yc;
public:

    FlowCylinder(
            const double D_,
            const double xc_,
            const double yc_,
            const double xMin_,
            const double xMax_,
            const double yMin_,
            const double yMax_,
            const double rho_,
            const double p0_,
            const double mu_,
            const double U0_
    ):FlowField(xMin_,xMax_,yMin_,yMax_,rho_,p0_,mu_,U0_),
      D(D_), xc(xc_), yc(yc_){}

    virtual double Psi(const double x_, const double y_) const;

    double ur(const double x_, const double y_) const;

    double uth(const double x_, const double y_) const;

    virtual double u(const double x_, const double y_) const;

    virtual double v(const double x_, const double y_) const;

    virtual Eigen::Matrix2d dVjdxi(const double x_, const double y_) const;

    virtual ~FlowCylinder() = default;
};



#endif //PARTICLETRAJECTORY_FLOWCYLINDER_H
