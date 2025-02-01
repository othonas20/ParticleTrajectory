#include "ParticleTrajectory.h"
//#include "VortexFlowField.h"
#include "FlowCylinder.h"
#include "FlowWall.h"
#include "Cylinder.h"
#include "WallStraight.h"
#include "FlowPipeLaminar.h"
#include "Pipe.h"
#include <Eigen/Dense>
#include <iostream>
#include <cstdlib> // for system()
#include <fstream>
#include <sstream>
#include <iomanip>  // for setprecision



#define sqr(x) pow(x,2)
#define PI 2 * acos(0)

int main()
{
    //// define parameters
        std::string Exercise = "B";
        std::string subExercise = "2a";

        const double rho = 1;
        const double rhop = 1000;
        const double mu = 1.8e-5;
        const double mup = 0.1;
        const double p0 = 0;

        const double U0 = 1;
        const double Rep = 0.5;
        const double Stk = 0.001;

        const double Dp = Rep * mup / (rhop * U0);
        const double D = 1./18. * Rep * mup/mu * Dp/Stk;
        const double Re = rho * U0 * D / mu;
        const double t0 = Stk * D / U0;

    //// construct solid wall
        const Eigen::Vector2d pointInternalWall {0,0};
        Cylinder myWall(pointInternalWall, 0.001*D, D);


    //// number of particles
        const int NumberParticles = 20;

        //// define initial conditions for every particle
        //const double A = U0 / (2*D);
        const double dx = 2 * D / (NumberParticles - 1);
        Eigen::Vector<Eigen::VectorXd, NumberParticles> initialCond;

        for (int i = 0; i < NumberParticles; ++i) {
            Eigen::VectorXd X0(6);

            X0(0) = - D;                         // x0
            X0(1) = - D + i * dx;                // y0
            X0(2) = 0;                           // u0
            X0(3) = 0;                           // v0
            X0(4) = (PI/6) * pow(Dp,3) * rhop;   // m0
            X0(5) = 300;                         // T0

            initialCond(i) = X0;
        }

    //// define solution parameters
        const double dt = t0 / 1; // = t0 / 20 for Runge Kutta, t0 / 1 for fix point
        const double tStop = 4000 * dt;
        const int intervalSaveData = int(tStop / dt) / 400; // = 5-10 for Runge Kutta, 1 for fix point

    //// particle counter
        int count = 0;

    //// construct particles and solving for its trajectories
        Eigen::Vector<ParticleTrajectory, NumberParticles> particles;

    //// Generate the result directory string
        // Use a stringstream to format with appropriate precision
        std::ostringstream oss;
        // Control the precision based on the value of Stk
        oss << std::fixed << std::setprecision(Stk == static_cast<int>(Stk) ? 0 : 1) << Stk;
        std::string resultDirectory = "Trajectories\\Exercise" + Exercise + "\\" + subExercise + "\\ST_" + oss.str();
    //// create the directory using mkdir
        std::string command = "mkdir " + resultDirectory + "\\data";
        // Execute the command
        int result = system(command.c_str());

        // Check if the command succeeded
        if (result == 0) {
            std::cout << "Directory created successfully: " << resultDirectory << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << resultDirectory << std::endl;
        }

        //// print parameters
        std::cout << "Re = " << Re << "\n";
        std::cout << "Dp = " << Dp << "\n";
        std::cout << "D = " << D << "\n";
        std::cout << "t0 = " << t0 << "\n";
        std::cout << "dt = " << dt << "\n";
        std::cout << "tStop = " << tStop << "\n";

        //// write parameters into file
        // Create and open the file
        std::ofstream outFile(resultDirectory + "\\parameters.txt");

        // Check if the file was created successfully
        if (!outFile) {
            std::cerr << "Error creating file 'parameters.txt'!" << std::endl;
        }

        // Write the parameters to the file
        outFile << "rho = " << rho << "\n";
        outFile << "rhop = " << rhop << "\n";
        outFile << "mu = " << mu << "\n";
        outFile << "mup = " << mup << "\n";
        outFile << "p0 = " << p0 << "\n";
        outFile << "U0 = " << U0 << "\n";
        outFile << "Rep = " << Rep << "\n";
        outFile << "Stk = " << Stk << "\n";
        outFile << "Dp = " << Dp << "\n";
        outFile << "D = " << D << "\n";
        outFile << "Re = " << Re << "\n";
        outFile << "t0 = " << t0 << "\n";

        // Close the file
        outFile.close();

        std::cout << "Parameters written to 'parameters.txt' successfully." << std::endl;

        //// construct flow field and store its data into a file
        //FlowPipeLaminar myField(D, 0, 5*D, -D/2, D/2, rho, p0, mu, U0);
        //FlowWall myField(A,-A/2,A/2,0,A,rho,p0,mu,U0);
        FlowCylinder myField(D,0,0,-4*D,2*D,-2*D,2*D,rho,p0,mu,U0);
        myField.PrintToFile(40, resultDirectory + "\\cylinderFlow_Stk_" + oss.str());


    //// run the ODE loops for every particle and save the results to files
        //for (int waveI = 0; waveI < 3; ++waveI) {
            for (int i = 0; i < NumberParticles; ++i) {
                std::cout<< count++ << "\n";

                // constructing particles
                    //const ParticleTrajectory dummy(&myField, false, false, initialCond(i), tStop, intervalSaveData, Dp, rhop, dt);
                    const ParticleTrajectory dummy(&myField,&myWall,true,false,initialCond(i),tStop,intervalSaveData,Dp,rhop,dt);
                    dummy.IncludeRandomMovement();
                    particles(i) = dummy;
                // solving with runge kutta loop
                    particles(i).LoopRungeKutta4th();
                    //particles(i).LoopImplicitEuler(1e-4);
                // storing results into files
                    particles(i).PrintToFile(resultDirectory +
                                         "//data//Trajectory_" + std::to_string(i));
            }
        //}

        //// create gnuplot file to visualise the flowfield and run it
        std::string pngName = "cylinderFlow";
        const double vectorLengthMultiplier = 0.5; // trial and error
        #include "plotFlowField.h"

            //// create gnuplot file to visualise the trajectories in a gif file and run it
        #include "plotTrajectory.h"


    return 0;
}
