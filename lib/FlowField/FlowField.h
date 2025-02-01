//
// Created by VS0121 on 18/1/2024.
//

#ifndef FLUIDVELOCITYFIELD_FLOWFIELDSTREAMFUNCTION_H
#define FLUIDVELOCITYFIELD_FLOWFIELDSTREAMFUNCTION_H

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "SolidWall.h"


// abstract class of the flow field
// the methods needed to be defined are:
//
//          virtual double Psi(const double x, const double y) const
//          virtual double u(const double x, const double y) const
//          virtual double v(const double x, const double y) const
//          virtual Eigen::Matrix2d dVjdxi(const double x, const double y)
class FlowField {
protected:

    // mutable cause they can change in function FillXYbounds
    mutable double xMin;
    mutable double xMax;
    mutable double yMin;
    mutable double yMax;
    double rho;
    double p0;
    double mu;
    double U0;

public:

    FlowField
            (
                    const double xMin_,
                    const double xMax_,
                    const double yMin_,
                    const double yMax_,
                    const double rho_,
                    const double p0_,
                    const double mu_,
                    const double U0_
            );

    virtual double Psi(const double x, const double y) const = 0;
    virtual double u(const double x, const double y) const = 0;
    virtual double v(const double x, const double y) const = 0;
    virtual Eigen::Matrix2d dVjdxi(const double x, const double y) const = 0;

    virtual double p(const double x, const double y) const;

    virtual void PrintToFile(const int N, const std::string filename) const;

    virtual void PrintToFile(const int N, const std::string filename, SolidWall* solidwall_pntr, const double distanceForTconst) const;

    void FillXYbounds(const double xmin, const double xmax, const double ymin, const double ymax) const;

    double getXmax() const;
    double getXmin() const;
    double getYmax() const;
    double getYmin() const;
    double getRho() const;
    double getP0() const;
    double getmu() const;
    double getU0() const;

    virtual ~FlowField() = default;

    double T_flowField(const double x_, const double y_, SolidWall *solidwall_pntr, const double distanceForTconst) const;
};





#endif //FLUIDVELOCITYFIELD_FLOWFIELDSTREAMFUNCTION_H
