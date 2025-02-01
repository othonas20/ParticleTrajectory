//
// Created by VS0121 on 19/1/2024.
//

#include "ParticleTrajectory.h"

// constructor without solid wall
ParticleTrajectory::ParticleTrajectory(FlowField* flowField_pntr_,
                                        const bool solid_wall_exist_,
                                        const bool temperature_and_humidity_change_,
                                        const Eigen::VectorXd& initialValuesDependent_,
                                        const double tStop_,
                                        const int dtSave_,
                                        const double D_,
                                        const double rho_p_,
                                        const double dt_)
                    :flowField_pntr(flowField_pntr_),
                    solidwall_pntr(nullptr),
                    solid_wall_exist(solid_wall_exist_),
                    include_random_movement(false),
                    distanceForTconst(0),
                    tStop(tStop_),
                    temperature_and_humidity_change(temperature_and_humidity_change_),
                    D(D_),
                    rho_p(rho_p_),
                    ODE(6,dt_, dtSave_, initialValuesDependent_, 0){}


// constructor with solid wall
ParticleTrajectory::ParticleTrajectory(FlowField* flowField_pntr_,
                                        SolidWall* solidwall_pntr_,
                                        const bool solid_wall_exist_,
                                        const bool temperature_and_humidity_change_,
                                        const Eigen::VectorXd& initialValuesDependent_,
                                        const double tStop_,
                                        const int dtSave_,
                                        const double D_,
                                        const double rho_p_,
                                        const double dt_)
                    :flowField_pntr(flowField_pntr_),
                    solidwall_pntr(solidwall_pntr_),
                    solid_wall_exist(solid_wall_exist_),
                    include_random_movement(false),
                    distanceForTconst(0),
                    tStop(tStop_),
                    temperature_and_humidity_change(temperature_and_humidity_change_),
                    D(D_),
                    rho_p(rho_p_),
                    ODE(6, dt_, dtSave_, initialValuesDependent_, 0){}

void ParticleTrajectory::CheckForDeflection(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const
{
    if (solid_wall_exist) {
        double x_ = dependentOldValues(0);
        double y_ = dependentOldValues(1);
        double F_ = solidwall_pntr->F(x_,y_);
        double F_crit = solidwall_pntr->Fcrit(x_,y_);
        double u_p_ = dependentOldValues(2);
        double v_p_ = dependentOldValues(3);

        if (F_ < F_crit)
        {
            Eigen::Vector2d normal = solidwall_pntr->Normal(x_, y_);
            Eigen::Vector2d tangent = solidwall_pntr->Tangent(x_, y_);

            Eigen::Vector2d V0(u_p_, v_p_);

            if (V0.dot(tangent) < 0)
                tangent *= -1;

            Eigen::Vector2d Vnt = - V0.dot(normal) * normal
                                  + V0.dot(tangent) * tangent;


            dependentOldValues(2) = Vnt(0);
            dependentOldValues(3) = Vnt(1);

//              dependentOldValues(2) *= -1;
//              dependentOldValues(3) *= -1;
        }
    }

}

void ParticleTrajectory::CheckForOutOfBounds(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const
{
    double x_ = dependentOldValues(0);
    double y_ = dependentOldValues(1);

    const double xMin = flowField_pntr->getXmin();
    const double xMax = flowField_pntr->getXmax();
    const double yMin = flowField_pntr->getYmin();
    const double yMax = flowField_pntr->getYmax();


    if (x_ < 0.99 * xMin) {
        if (Stop_Loop_If_Out_Of_Bounds) independentOldValue = 2 * tStop;
        else dependentOldValues(0) = 0.98 * xMax;
    }
    else if (x_ > 0.99 * xMax) {
        if (Stop_Loop_If_Out_Of_Bounds) independentOldValue = 2 * tStop;
        else dependentOldValues(0) = 0.98 * xMin;
    }
    else if (y_ < 0.99 * yMin) {
        if (Stop_Loop_If_Out_Of_Bounds) independentOldValue = 2 * tStop;
        else dependentOldValues(1) = 0.98 * yMax;
    }
    else if (y_ > 0.99 * yMax) {
        if (Stop_Loop_If_Out_Of_Bounds) independentOldValue = 2 * tStop;
        else dependentOldValues(1) = 0.98 * yMin;
    }

}

void ParticleTrajectory::PrintToFile(const std::string filename) const
{
    std::ofstream file(filename, std::ios::trunc);
    if (file.is_open())
    {
        std::cout << "file opened!\n";
        file << "#t   x   y   z   u   v   m   T\n";

        const Eigen::RowVectorXd& t = varIndependent;
        const Eigen::RowVectorXd& x = varDependent.row(0);
        const Eigen::RowVectorXd& y = varDependent.row(1);
        const Eigen::RowVectorXd& u = varDependent.row(2);
        const Eigen::RowVectorXd& v = varDependent.row(3);
        const Eigen::RowVectorXd& m = varDependent.row(4);
        const Eigen::RowVectorXd& T = varDependent.row(5);

        const int N = t.size();
        for (int i = 0; i < N ; ++i) {
            file << t(i) << " " << x(i) << " " << y(i) << " " << 0 << " " << u(i) << " " << v(i) << " " << m(i) << " " << T(i) <<"\n";
        }
        file.close();
    }
    else
    {
        std::cout << " file not opened!\n";
    }
}

const FlowField* ParticleTrajectory::FlowFieldPtr() const
{
    return flowField_pntr;
}

// Pure virtual functions inherited from ODE class

bool ParticleTrajectory::TerminatingCondition(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const
{
    CheckForDeflection(independentOldValue, dependentOldValues);
    CheckForOutOfBounds(independentOldValue, dependentOldValues);
    return independentOldValue < tStop && dependentOldValues(4) / this->varDependent(4,0) > 1e-2;
}

const Eigen::VectorXd ParticleTrajectory::RightSideODEs(const double independentOldValue, const Eigen::VectorXd& dependentOldValues) const
{
    const double PI = 3.14159265358979323846;
    const double x_ = dependentOldValues(0);
    const double y_ = dependentOldValues(1);
    double u_p_ = dependentOldValues(2);
    double v_p_ = dependentOldValues(3);

    if (include_random_movement)
    {
        const double Vmean = sqrt(u_p_*u_p_ + v_p_*v_p_);
        const double Vstddev = 0.05 * Vmean;

        // Create a normal (Gaussian) distribution with the given mean and standard deviation
        std::normal_distribution<> distribution_up(u_p_, Vstddev);
        std::normal_distribution<> distribution_vp(v_p_, Vstddev);

        // Create a random number generator
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

        u_p_ = distribution_up(gen);
        v_p_ = distribution_vp(gen);
    }

    const double m_p_ = dependentOldValues(4);
    const double T_p_ = dependentOldValues(5);
    const double Vol_p_ = m_p_ / rho_p;
    const double D_p_ = pow(6 * Vol_p_ / PI, 0.3333333333);
    const double m_f_ = flowField_pntr->getRho() * Vol_p_;
    const double mu_f = flowField_pntr->getmu();
    const Eigen::Vector2d g(0,-9.81);
    const Eigen::Vector2d V_p_(u_p_, v_p_);
    const Eigen::Vector2d V_f_(flowField_pntr->u(x_,y_), flowField_pntr->v(x_,y_));
    const Eigen::Matrix2d gradV_f(flowField_pntr->dVjdxi(x_,y_));

    const Eigen::Vector2d V_p_dot_grad_V_f = {V_p_.dot(gradV_f.col(0)),
                                               V_p_.dot(gradV_f.col(1))};
    const Eigen::Vector2d V_f_dot_grad_V_f = {V_f_.dot(gradV_f.col(0)),
                                              V_f_.dot(gradV_f.col(1))};

    const double Rer = flowField_pntr->getRho() * (V_p_ - V_f_).norm() * D / flowField_pntr->getmu();
    double CD = (Rer < 800) ? ((24 / Rer) * (1 + 0.15 * pow(Rer, 0.687)) ) : 0.45;
    if (Rer < 1) { CD = 0.45 / Rer;}
    Eigen::VectorXd returnVec(6);
    const Eigen::Vector2d dV_p_dt = (
                                            (m_p_ - m_f_) * g // gravity + elevation
                                            + 0.5 * CD * flowField_pntr->getRho() * 0.25 * PI * pow(D,2) * (V_p_ - V_f_).norm() * (V_f_ - V_p_) // aerodynamic drag
                                            + 0.5 * m_f_ * V_p_dot_grad_V_f // additional mass
                                            + m_f_ * V_f_dot_grad_V_f // undisturbed flow
                                    ) / (m_p_ + 0.5 * m_f_); // devide with total mass

    returnVec(0) = u_p_; // dx/dt = up
    returnVec(1) = v_p_; // dx/dt = vp
    returnVec(2) = dV_p_dt.x(); // return the x value of the resulted vector
    returnVec(3) = dV_p_dt.y(); // return the y value of the resulted vector

    if (!temperature_and_humidity_change){
        returnVec(4) = 0;
        returnVec(5) = 0;
    }else
    {
        const double Sh = 2;
        const double Dv = 5e-9;
        const double RH_inf = 0.50;
        const double r = 0.005;
        const double Mc = 28.89;
        const double Md = 18;
        const double M = (Mc + r*Md) / (1+r);
        const double A = 8.07131;
        const double B = 1730.63;
        const double C = 233.426;
        auto Psat = [A, B, C](double T) -> double { return 133.322 * std::pow(10.0, A - (B / (T + C)));};
        const double P = FlowFieldPtr()->p(x_,y_);
        const double Pinf = FlowFieldPtr()->getP0();
        const double Tinf = 300;
        const double Psat_inf = Psat(Tinf - 273.15);
        const double omega_s_inf = (Md / M) * (Psat_inf / Pinf);
        const double omega_inf = RH_inf * omega_s_inf;
        const double Psat_ = Psat(T_p_ - 273.15);
        const double omega_s = (Md / M) * (Psat_ / P);

        const double Cp = 4186;
        const double k = 0.598;
        const double mu = 0.001;
        const double Pr = Cp * mu / k;
        const double Nu = 2 + 0.4 * pow(Rer, 0.5) * pow(Pr, 0.3333);
        const double T_flowfield = flowField_pntr->T_flowField(x_,y_, solidwall_pntr, distanceForTconst);
        const double Area = 0.25 * PI * pow(D_p_,2);
        const double h = Nu * k / D;
        const double Qin = h * Area * (T_flowfield - T_p_);
        const double Lamda = 2260000; // J/kg
        double dmdt = (T_p_ < 373) ? Sh * flowField_pntr->getRho() * PI * D_p_ * Dv * (omega_inf - omega_s) : - Qin/Lamda;
        double dTdt = (T_p_ < 373) ? (Qin - dmdt * Cp * T_p_) / (m_p_ * Cp) : 0;
        returnVec(4) = dmdt; // dm/dt
        returnVec(5) = dTdt; // dT/dt
    }

    return returnVec;
}


