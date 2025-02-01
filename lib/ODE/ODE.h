//
// Created by VS0121 on 4/6/2023.
//

#ifndef ODE_H
#define ODE_H

#include <iostream>
#include <Eigen/Dense>


// Abstract class, to solve a system of <numberOfEquations> ODE's
class ODE
{

protected:

// Member variables

    int numberOfEquations;

    // Matrix storing the dependent variables
    mutable Eigen::MatrixXd varDependent;

    // Vector storing the independent variable
    mutable Eigen::RowVectorXd varIndependent;

    // Step size of the integration
    double interval;

    // save interva;
    int tSave;

public:


    // Constructor
    ODE
    (
            int numberOfEquations_,
            double interval_,
            const int tSave_,
            const Eigen::VectorXd& initialValuesDependent_,
            double initialValueIndependent_
    );

    // default constructor
    ODE() = default;

    // Pure virtual function to be implemented by derived classes
    // It should fill the right-hand side of the ODEs for the specific problem
    virtual const Eigen::VectorXd RightSideODEs(const double independentOldValue, const Eigen::VectorXd& dependentOldValues) const = 0;

    virtual void LoopImplicitEuler(const double err) const;

    // Method to perform a solving-loop, using the Runge-Kutta 4th order method
    virtual void LoopRungeKutta4th() const;

    // Pure virtual function for the terminating condition to be implemented by derived classes
    // It should determine when to terminate the integration
    virtual bool TerminatingCondition(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const = 0;

    const Eigen::RowVectorXd& getIndependentVariable() const;

    const Eigen::MatrixXd& getDependentVariable() const;

    // defualt destructor
    virtual ~ODE() = default;


};




#endif //ODE_H
