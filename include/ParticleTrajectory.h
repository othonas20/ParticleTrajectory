//
// Created by VS0121 on 12/12/2023.
//

#ifndef PARTICLETRAJECTORY_PARTICLETRAJECTORY_H
#define PARTICLETRAJECTORY_PARTICLETRAJECTORY_H

#include <memory>
#include <fstream>
#include <random>
#include "FlowField.h"
//#include "VortexFlowField.h"
#include "SolidWall.h"
#include "ODE.h"



class ParticleTrajectory : public ODE{

private:

//    std::unique_ptr<FlowField> flowField_pntr;
//    std::unique_ptr<SolidWall> solidwall_pntr;
    FlowField* flowField_pntr;
    SolidWall* solidwall_pntr;
    mutable double D; // particle diametre in m
    mutable double rho_p; // particle density
    bool solid_wall_exist;
    bool temperature_and_humidity_change;
    mutable bool include_random_movement;
    mutable bool Stop_Loop_If_Out_Of_Bounds = false;
    mutable double distanceForTconst;
    double tStop;

public:

    // constructor without solid wall
    // ParticleTrajectory(FlowField* flowField_pntr_,
    //    bool solid_wall_exist_,
    //    bool temperature_and_humidity_change_,
    //    const Eigen::VectorXd& initialValuesDependent_,  !!! MUST BE OF SIZE 6 !!!
    //    const double tStop_,
    //    const int dtSave_,
    //    const double D_,
    //    const double rho_p_,
    //    const double dt_);
    ParticleTrajectory(FlowField* flowField_pntr_,
                       bool solid_wall_exist_,
                       bool temperature_and_humidity_change_,
                           const Eigen::VectorXd& initialValuesDependent_,
                       const double tStop_,
                       const int dtSave_,
                       const double D_,
                       const double rho_p_,
                       const double dt_);

    // // constructor with solid wall
    // ParticleTrajectory(FlowField* flowField_pntr_,
    //                    SolidWall* solidwall_pntr_,
    //                    const bool solid_wall_exist_,
    //                    const bool temperature_and_humidity_change_,
    //                    const Eigen::VectorXd& initialValuesDependent_,    !!! MUST BE OF SIZE 6 !!!
    //                    const double tStop_,
    //                    const int dtSave_,
    //                    const double D_,
    //                    const double rho_p_,
    //                    const double dt_);
    ParticleTrajectory(FlowField* flowField_pntr_,
                       SolidWall* solidwall_pntr_,
                       const bool solid_wall_exist_,
                       const bool temperature_and_humidity_change_,
                       const Eigen::VectorXd& initialValuesDependent_,
                       const double tStop_,
                       const int dtSave_,
                       const double D_,
                       const double rho_p_,
                       const double dt_);

    // default constructor
    ParticleTrajectory() = default;

    void CheckForDeflection(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const;

    void CheckForOutOfBounds(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const;

    void IncludeRandomMovement() const {include_random_movement = true;};

    void StopLoopIfOutOfBounds() const {Stop_Loop_If_Out_Of_Bounds = true;};

    void SetDistanceForTconst(const double d) const {distanceForTconst = d;}

    void PrintToFile(const std::string filename) const;

    const FlowField* FlowFieldPtr() const;

// Pure virtual functions inherited from ODE class

    virtual bool TerminatingCondition(double& independentOldValue, Eigen::VectorXd& dependentOldValues) const;

    virtual const Eigen::VectorXd RightSideODEs(const double independentOldValue, const Eigen::VectorXd& dependentOldValues) const;

};

double T_flowField(const double x_, const double y_);

#endif //PARTICLETRAJECTORY_PARTICLETRAJECTORY_H
