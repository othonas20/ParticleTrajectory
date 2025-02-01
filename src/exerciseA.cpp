//
// Created by VS0121 on 22/8/2024.
//
#include "ParticleTrajectory.h"
#include "VortexFlowField.h"
#include "FlowCylinder.h"
#include "Cylinder.h"
#include <Eigen/Dense>
#include <iostream>
#include <cstdlib> // for system()
#include <fstream>
#include <sstream>
#include <iomanip>  // for setprecision

#define sqr(x) pow(x,2)
#define PI 2 * acos(0)


int main() {

    //// define parameters
    std::string Exercise = "A";
    std::string subExercise = "2";

    const double rho = 1;
    const double mu = 1.8e-5;
    const double p0 = 0;

    const double R = 1e-6;
    const double W = 0.5;
    const double Stk = 0.001;
    const double A = 1 / Stk;
    const double Rep = .1;

    //// Generate the result directory string
    // Use a stringstream to format with appropriate precision
    std::ostringstream oss;
    // Control the precision based on the value of Stk
    oss << std::fixed << std::setprecision(Stk == static_cast<int>(Stk) ? 0 : 1) << Stk;
    std::string resultDirectory = "Trajectories/Exercise" + Exercise + "/" + subExercise + "/ST_" + oss.str();
    //std::string resultDirectory = "Trajectories\\Exercise" + Exercise + "\\" + subExercise + "\\ST_" + oss.str();
    
    //// create the directory using mkdir
    std::string command = "mkdir -p " + resultDirectory + "/data";
    //std::string command = "mkdir " + resultDirectory + "\\data";
    
    // Execute the command
    int result = system(command.c_str());

    // Check if the command succeeded
    if (result == 0) {
        std::cout << "Directory created successfully: " << resultDirectory << std::endl;
    } else {
        std::cerr << "Failed to create directory: " << resultDirectory << std::endl;
    }

    //// construct flow field and store its data into a file
    VortexFlowField myField(rho,p0,mu,R,Rep,W,A);

    const double a = myField.get_a();
    Eigen::Vector2d center {0.5 * PI * a, PI * a};
    double xMin = -1 * PI * a + center.x();
    double xMax = 1 * PI * a + center.x();
    double yMin = -1 * PI * a + center.y();
    double yMax = 1 * PI * a + center.y();
    myField.FillXYbounds(xMin,xMax,yMin,yMax);
    myField.PrintToFile(40, resultDirectory + "/vortexFlow_Stk_" + oss.str());
    //myField.PrintToFile(50, resultDirectory + "\\vortexFlow_Stk_" + oss.str());

    //// particle parameters
    const int NumberParticles = 100;
    const double D = Rep * mu / (rho * myField.getU0());
    const double rhop = (1/R - 0.5) * rho;

    //// define initial conditions for every particle
    const double dx = (xMax - xMin) / (NumberParticles - 1);
    Eigen::Vector<Eigen::VectorXd, NumberParticles> initialCond;

    for (int i = 0; i < NumberParticles; ++i) {
        Eigen::VectorXd X0(6);

        X0(0) = xMin + i * dx;                      // x0
        X0(1) = 0.98 * yMax;               // y0
        X0(2) = 0;                           // u0
        X0(3) = 0;                           // v0
        X0(4) = (PI/6) * pow(D,3) * rhop;   // m0
        X0(5) = 300;                         // T0

        initialCond(i) = X0;
    }

    //// define solution parameters
    const double t0 = rhop * pow(D, 2) / (18 * mu);

    const double dt = t0 / 2; // = t0 / 20 for Runge Kutta, t0 / 1 for fix point
    const int Nstop = 500000;
    const int Nsave = 200;
    const double tStop = Nstop * dt;
    const int intervalSaveData = int(Nstop / Nsave);

    //// print parameters
    std::cout << "dt = " << dt << "\n";
    std::cout << "tStop = " << tStop << "\n";

    //// write parameters into file
    // Create and open the file
    std::ofstream outFile(resultDirectory + "/parameters.txt");
    //std::ofstream outFile(resultDirectory + "\\parameters.txt");

    // Check if the file was created successfully
    if (!outFile) {
        std::cerr << "Error creating file 'parameters.txt'!" << std::endl;
    }

    // Write the parameters to the file
    outFile << "rho = " << rho << "\n";
    outFile << "rhop = " << rhop << "\n";
    outFile << "mu = " << mu << "\n";
    outFile << "p0 = " << p0 << "\n";
    outFile << "U0 = " << myField.getU0() << "\n";
    outFile << "Rep = " << Rep << "\n";
    outFile << "Stk = " << Stk << "\n";
    outFile << "Dp = " << D << "\n";
    outFile << "Re = " << rho * myField.getU0() * myField.get_a() / mu<< "\n";
    outFile << "t0 = " << t0 << "\n";

    // Close the file
    outFile.close();

    std::cout << "Parameters written to 'parameters.txt' successfully." << std::endl;


    //// particle counter
    int count = 0;

    //// construct particles and solving for its trajectories
    Eigen::Vector<ParticleTrajectory, NumberParticles> particles;

    // run the ODE loops for every particle and save the results to files
    for (int i = 0; i < NumberParticles; ++i) {
        std::cout<< count++ << "\n";

        // constructing particles
        const ParticleTrajectory dummy(&myField, false, false, initialCond(i), tStop, intervalSaveData, D, rhop, dt);
        particles(i) = dummy;
        // solving with runge kutta loop
        particles(i).LoopRungeKutta4th();
        // storing results into files
        particles(i).PrintToFile(resultDirectory +
                                 "/data/Trajectory_" + std::to_string(i));
                                 //"//data//Trajectory_" + std::to_string(i));
    }

    //// create gnuplot file to visualise the flowfield and run it
    const std::string pngName = "vortexFlow";
    const double vectorLengthMultiplier = 5; // trial and error
    #include "plotFlowField.h"

    //// create gnuplot file to visualise the trajectories in a gif file and run it
    #include "plotTrajectory.h"



    return 0;
}
