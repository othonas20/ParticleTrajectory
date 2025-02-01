//
// Created by VS0121 on 19/8/2024.
//

#include "FlowCylinder.h"

    double FlowCylinder::Psi(const double x_, const double y_) const
    {
        const double r = sqrt(pow((x_-xc),2) + pow((y_-yc),2));
        const double theta = atan2(y_-yc, x_-xc);
        return U0 * r * sin(theta) * (1 - 0.25 * pow(D/r, 2));
    }

    double FlowCylinder::ur(const double x_, const double y_) const
    {
        const double r = sqrt(pow((x_-xc),2) + pow((y_-yc),2));
        const double theta = atan2(y_-yc, x_-xc);
        return (r > D/2) ? U0 * cos(theta) * (1 - 0.25 * pow(D/r,2)) : 0;
    }

    double FlowCylinder::uth(const double x_, const double y_) const
    {
        const double r = sqrt(pow((x_-xc),2) + pow((y_-yc),2));
        const double theta = atan2(y_-yc, x_-xc);
        return (r > D/2) ?  - U0 * sin(theta) * (1 + 0.25 * pow(D/r, 2)) : 0;
    }

    double FlowCylinder::u(const double x_, const double y_) const
    {
        const double theta = atan2(y_ - yc, x_ - xc);
        return ur(x_, y_) * cos(theta) - uth(x_,y_) * sin(theta);
    }

    double FlowCylinder::v(const double x_, const double y_) const
    {
        const double theta = atan2(y_ - yc, x_ - xc);
        return ur(x_, y_) * sin(theta) + uth(x_,y_) * cos(theta);
    }

    Eigen::Matrix2d FlowCylinder::dVjdxi(const double x_, const double y_) const
    {
        const double e = 1e-4;
        Eigen::Matrix2d gradV;
        gradV(0,0) = (u(x_ + e, y_) - u(x_ - e, y_)) / (2*e); // du/dx
        gradV(0,1) = (v(x_ + e, y_) - v(x_ - e, y_)) / (2*e); // dv//dx
        gradV(1,0) = (u(x_, y_ + e) - u(x_, y_ - e))/ (2*e); // du/dy
        gradV(1,1) = (v(x_, y_ + e) - v(x_, y_ - e)) / (2*e); // dv/dy

        return gradV;
    }



