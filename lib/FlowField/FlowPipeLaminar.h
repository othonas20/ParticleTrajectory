//
// Created by VS0121 on 23/8/2024.
//

#ifndef PARTICLETRAJECTORY_FLOWPIPELAMINAR_H
#define PARTICLETRAJECTORY_FLOWPIPELAMINAR_H

#include "FlowField.h"

class FlowPipeLaminar : public FlowField{
private:
    double D;
    double Umax;

public:

    FlowPipeLaminar(const double D_,
                    const double xMin_,
                    const double xMax_,
                    const double yMin_,
                    const double yMax_,
                    const double rho_,
                    const double p0_,
                    const double mu_,
                    const double U0_
                    ): FlowField(xMin_,xMax_,yMin_,yMax_,rho_,p0_,mu_,U0_),
                    D(D_), Umax(2 * U0_){}

    virtual double Psi(const double x_, const double y_) const {return 0;};

    virtual double u(const double x_, const double y_) const;

    virtual double v(const double x_, const double y_) const;

    virtual Eigen::Matrix2d dVjdxi(const double x_, const double y_) const;

    virtual double p(const double x_, const double y_) const override;

    double getUmax() const { return Umax;};

    bool PointInsideBounds(const double x_, const double y_) const;

    virtual ~FlowPipeLaminar() = default;



};


#endif //PARTICLETRAJECTORY_FLOWPIPELAMINAR_H
