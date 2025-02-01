#pragma once
#include "FlowField.h"

class VortexFlowField : public FlowField{
private:
    double a;
public:
    VortexFlowField(
        const double a_,
        const double U0_, 
        const double xMin_, 
        const double xMax_, 
        const double yMin_, 
        const double yMax_, 
        const double rho_, 
        const double p0_, 
        const double mu_
        )
                        : FlowField
                        (
                        xMin_,
                        xMax_,
                        yMin_,
                        yMax_,
                        rho_,
                        p0_,
                        mu_,
                        U0_
                        ),
                        a(a_){}

    VortexFlowField(const double rho_,
                    const double p0_,
                    const double mu_,
                    const double R,
                    const double Re,
                    const double W,
                    const double A
        );


    virtual double Psi(const double x, const double y) const;

    virtual double u(const double x, const double y) const;
    
    virtual double v(const double x, const double y) const;

    virtual Eigen::Matrix2d dVjdxi(const double x, const double y) const;

    double get_a() const {return a;};

    virtual ~VortexFlowField() = default;
};

