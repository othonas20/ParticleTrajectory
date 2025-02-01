#include "VortexFlowField.h"

    VortexFlowField::VortexFlowField(const double rho_, const double p0_, const double mu_, const double R, const double Re,
                                     const double W, const double A) : FlowField(0,0,0,0,rho_,p0_,mu_, 0), a(0){
        std::cout << "Dont forget to initialise x-y min-max !\n";
        const double rhop = (1/R - 0.5) * rho;
        U0 = pow(
                abs((1./R) - 1.5) * mu * pow(Re,2) / (3 * rho * W)
                ,
                0.33333333333
        );
        const double D = Re * mu / (rho * U0);
        a = A * U0 * pow(D,2) * (rhop + 0.5 * rho) / (18 * mu);

        const double Stk = 1 / A;
        const double t0 = rhop * powf(D, 2) / (18 * mu);
        const double V_p = 3.14 * powf(D, 3) / 6;
        const double mp = rhop * V_p;

        const double mf = rho * V_p;
        const double Ws = (mp - mf) * 9.81 / (3 * 3.14 * D * mu);
        std::cout << "r/a = " << 0.5 * D / a << "\n";
        std::cout << "r * Ws / nu = " << 0.5 * D * Ws / (mu / rho) << "\n";
        std::cout << "r^2 * U0 / (a * nu) = " << powf(D/2, 2) * U0 / (a * (mu / rho)) << "\n\n";
        std::cout << "a = " << a << "\n";


        // print particle parameters
        std::cout << "R = " << R << "\n";
        std::cout << "A = " << A << "\n";
        std::cout << "W = " << W << "\n";
        std::cout<< "rhop = " << rhop << "\n";
        std::cout<< "Stk = " << Stk << "\n";
        std::cout<< "Re particle = " << Re << "\n";
        std::cout<< "Re flow field = " << rho * U0 * a / mu << "\n";
        std::cout << "U0 = " << U0 << "\n";
        std::cout<< "D = " << D << "\n";
        std::cout<< "t0 = " << t0 << "\n";
        std::cout<< "V_p = " << V_p << "\n";
        std::cout<< "mp = " << mp << "\n";
    }

    double VortexFlowField::Psi(const double x, const double y) const
    {
        return a * U0 * cos(x / a) * sin(y / a);
    }

    double VortexFlowField::u(const double x, const double y) const
    {
        return U0 * cos(x / a) * cos(y / a);
    }

    double VortexFlowField::v(const double x, const double y) const
    {
        return U0 * sin(x / a) * sin(y / a);
    }

    Eigen::Matrix2d VortexFlowField::dVjdxi(const double x, const double y) const
    {
        Eigen::Matrix2d gradV;
        gradV(0,0) = - (U0 / a) * sin(x / a) * cos(y / a);
        gradV(0,1) =   (U0 / a) * cos(x / a) * sin(y / a);
        gradV(1,0) = - (U0 / a) * cos(x / a) * sin(y / a);
        gradV(1,1) =   (U0 / a) * sin(x / a) * cos(y / a);

        return gradV;
    }


