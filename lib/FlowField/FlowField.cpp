//
// Created by VS0121 on 18/1/2024.
//

#include "FlowField.h"
#include <fstream>

FlowField::FlowField
        (
                const double xMin_,
                const double xMax_,
                const double yMin_,
                const double yMax_,
                const double rho_,
                const double p0_,
                const double mu_,
                const double U0_
        )
        :xMin(xMin_),
        xMax(xMax_),
        yMin(yMin_),
        yMax(yMax_),
        rho(rho_),
        p0(p0_),
        mu(mu_),
        U0(U0_){}

double FlowField::p(const double x, const double y) const
{
    return p0 + 0.5 * rho * powf(U0,2) - 0.5 * rho * (powf(u(x,y),2) + powf(v(x,y),2));
}

void FlowField::PrintToFile(const int N, const std::string filename) const
{
    double dx = (xMax - xMin) / (N - 1);
    double dy = (yMax - yMin) / (N - 1);
    double X, Y, PSI, U, V, P;

    std::ofstream file(filename, std::ios::trunc);
    if (file.is_open())
    {
        std::cout << "file opened!\n";
        file << "x y z psi u v p\n";
        for (int i = 0; i < N ; ++i) {
            for (int j = 0; j < N ; ++j) {
                X = i * dx + xMin;
                Y = j * dy + yMin;
                PSI = Psi(X,Y);
                U = u(X,Y);
                V = v(X,Y);
                P = p(X,Y);

                file <<X << " " << Y << " " << 0 << " " << PSI << " " << U << " " << V << " " << P <<"\n";
            }
        }
        file.close();
    }
    else
    {
        std::cout << " file not opened!\n";
    }

}


void FlowField::PrintToFile(const int N, const std::string filename, SolidWall* solidwall_pntr, const double distanceForTconst) const {
    double dx = (xMax - xMin) / (N - 1);
    double dy = (yMax - yMin) / (N - 1);
    double X, Y, PSI, U, V, P, T;

    std::ofstream file(filename, std::ios::trunc);
    if (file.is_open())
    {
        std::cout << "file opened!\n";
        file << "x y z psi u v p T\n";
        for (int i = 0; i < N ; ++i) {
            for (int j = 0; j < N ; ++j) {
                X = i * dx + xMin;
                Y = j * dy + yMin;
                PSI = Psi(X,Y);
                U = u(X,Y);
                V = v(X,Y);
                P = p(X,Y);
                T = T_flowField(X, Y, solidwall_pntr, distanceForTconst);

                file <<X << " " << Y << " " << 0 << " " << PSI << " " << U << " " << V << " " << P << " " << T  <<"\n";
            }
        }
        file.close();
    }
    else
    {
        std::cout << " file not opened!\n";
    }
}

void FlowField::FillXYbounds(const double xmin, const double xmax, const double ymin,
                                           const double ymax) const
{
    xMin = xmin;
    xMax = xmax;
    yMin = ymin;
    yMax = ymax;
}

double FlowField::T_flowField(const double x_,
                              const double y_,
                              SolidWall* solidwall_pntr,
                              const double distanceForTconst) const
{
    const double Thigh = 0.75 * 900;
    const double Tlow = 300;
    const double d = solidwall_pntr->DistanceFromWall(x_,y_);
    if (solidwall_pntr->PointInFluidRegion(x_,y_))
    {
        return (d < distanceForTconst) ?
               Thigh - (Thigh - Tlow) * d / distanceForTconst
                                       : Tlow;
    }
    else return 900;

}

// Get the maximum value of x
double FlowField::getXmax() const
{
    return xMax;
}

// Get the minimum value of x
double FlowField::getXmin() const
{
    return xMin;
}

// Get the maximum value of y
double FlowField::getYmax() const
{
    return yMax;
}

// Get the minimum value of y
double FlowField::getYmin() const
{
    return yMin;
}

// Get the value of density (rho)
double FlowField::getRho() const
{
    return rho;
}

// Get the reference pressure (p0)
double FlowField::getP0() const
{
    return p0;
}

// Get the dynamic viscosity (mu)
double FlowField::getmu() const
{
    return mu;
}

// Get the characteristic velocity (U0)
double FlowField::getU0() const
{
    return U0;
}

