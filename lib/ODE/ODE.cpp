//
// Created by VS0121 on 4/6/2023.
//

#include "ODE.h"
#include <chrono>


// Constructor definition
ODE::ODE
        (
                const int numberOfEquations_,
                const double interval_,
                const int tSave_,
                const Eigen::VectorXd& initialValuesDependent_,
                const double initialValueIndependent_
        ):
        numberOfEquations(numberOfEquations_),
        interval(interval_),
        tSave(tSave_)
{
    // initialization of independent variable
    varIndependent.resize(1);
    varIndependent(0) = initialValueIndependent_;

    // initialization of dependent variables matrix
    varDependent = initialValuesDependent_;
    std::cout << "Runge Kutta 4th object created!\n";
}

void ODE::LoopImplicitEuler(const double err) const
{
    static int numberOfIterationsMean = 0;

    int SIZE = 1; // Current size of the varDependent matrix
    int t = 0;

    Eigen::VectorXd oldValuesDependent(numberOfEquations); // Values of the dependent variables at the current step
    oldValuesDependent = varDependent.col(SIZE - 1);
    double oldValueIndependent = varIndependent(SIZE - 1); // Value of the independent variable at the current step

    std::cout<<"Beggining Fix Point algorithm loop:\n";

    while (TerminatingCondition(oldValueIndependent, oldValuesDependent))
    {

        double residual = 1;
        Eigen::VectorXd dY; // = dx * f(told, y)

        // update old values for dependent variables
        while (residual > err)
        {
            numberOfIterationsMean++;

            dY = interval * RightSideODEs(oldValueIndependent, oldValuesDependent);
            residual = dY.norm() / oldValuesDependent.norm();
            oldValuesDependent += dY;
        }

        oldValueIndependent += interval; // update old value for independent variable

        t++;
        if (t == tSave){
            t = 0;
            SIZE++;

            // Resize and update varDependent matrix to store the updated values
            varDependent.conservativeResize(numberOfEquations, SIZE);
            varDependent.col(SIZE - 1) = oldValuesDependent;
            // Resize and update varIndependent vector for the next step
            varIndependent.conservativeResize(SIZE);
            varIndependent(SIZE - 1) = oldValueIndependent;
        }
    }
    std::cout<<"End of Runge - Kutta 4th order loop\n";
}


void ODE::LoopRungeKutta4th() const
{
    Eigen::MatrixXd K(numberOfEquations, 4);     // Matrix to store the intermediate values K1, K2, K3, K4

    int SIZE = 1; // Current size of the varDependent matrix
    int t = 0;

    Eigen::VectorXd oldValuesDependent(numberOfEquations); // Values of the dependent variables at the current step
    oldValuesDependent = varDependent.col(SIZE - 1);
    double oldValueIndependent = varIndependent(SIZE - 1); // Value of the independent variable at the current step

    std::cout<<"Beggining Runge - Kutta 4th order loop:\n";

    while (TerminatingCondition(oldValueIndependent, oldValuesDependent))
    {
        // Calculate the intermediate values K1, K2, K3, K4 using the right-hand side of the ODEs
        K.col(0) = RightSideODEs(oldValueIndependent,oldValuesDependent);
        K.col(1) = RightSideODEs(oldValueIndependent + interval/2  ,oldValuesDependent + 0.5 * interval * K.col(0));
        K.col(2) = RightSideODEs(oldValueIndependent + interval/2  ,oldValuesDependent + 0.5 * interval * K.col(1));
        K.col(3) = RightSideODEs(oldValueIndependent + interval    ,oldValuesDependent + 0.5 * interval * K.col(2));

        oldValuesDependent += interval * (K.col(0) + 2 * K.col(1) + 2 * K.col(2) + K.col(3)) / 6; // update old values for dependent variables
        oldValueIndependent += interval; // update old value for independent variable

        //std::cout << oldValueIndependent << "   " << oldValuesDependent.transpose() << std::endl;

        t++;
        if (t == tSave){
            t = 0;
            SIZE++;

            // Resize and update varDependent matrix to store the updated values
            varDependent.conservativeResize(numberOfEquations, SIZE);
            varDependent.col(SIZE - 1) = oldValuesDependent;
            // Resize and update varIndependent vector for the next step
            varIndependent.conservativeResize(SIZE);
            varIndependent(SIZE - 1) = oldValueIndependent;
        }
    }
    std::cout<<"End of Runge - Kutta 4th order loop\n";
}



const Eigen::RowVectorXd& ODE::getIndependentVariable() const
{
    return  varIndependent;
}

const Eigen::MatrixXd& ODE::getDependentVariable() const
{
    return varDependent;
}

