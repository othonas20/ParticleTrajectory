//
// Created by VS0121 on 22/8/2024.
//

// Create and open the plotTrajectory.p file
std::ofstream plotFile("plotTrajectory.p");

// Check if the file was created successfully
if (!plotFile) {
std::cerr << "Error creating plotTrajectory.p file!" << std::endl;
return 1;
}

// Write the provided gnuplot script into the file
plotFile << "# Set up the terminal and output file (GIF)\n"
<< "set terminal gif animate optimize delay 10 size 2000,2000\n"
<< "set output 'Trajectories/Exercise" + Exercise + "/" + subExercise + "/ST_" + oss.str() + "/trajectory_animation_C_ST_" + oss.str() + ".gif'\n\n"
<< "# Set plot properties\n"
<< "set title \"Trajectory Stk = " + oss.str() + "\"\n"
<< "set xlabel \"X\"\n"
<< "set ylabel \"Y\"\n"
<< "set grid\n\n"
<< "set size ratio -1\n"
<< "set colorbox\n\n"
<< "set pointsize 1\n"
<< "set style line 2 lc rgb '#0060ad' pt 7\n\n"
<< "# Define the animation parameters\n"
<< "file_path_flowfield = 'Trajectories/Exercise" + Exercise
<< + "/" + subExercise
<< + "/ST_" + oss.str() + "/" + pngName + "_Stk_" + oss.str() + "'\n\n"
<< "# Calculate the maximum magnitudes\n"
<< "stats file_path_flowfield using (sqrt($5**2 + $6**2)) name \"Velocities\"\n"
<< "stats file_path_flowfield using 1 name \"X\"\n"
<< "stats file_path_flowfield using 2 name \"Y\"\n\n"
<< "max_velocity = Velocities_max\n"
<< "DX = X_max - X_min\n"
<< "DY = Y_max - Y_min\n\n"
<< "set xrange [X_min - 0.2*DX : X_max + 0.2 * DX]\n"
<< "set yrange [Y_min - 0.2 * DY : Y_max + 0.2 * DY]\n\n"
<< "# Set the color palette range based on the maximum magnitude\n"
<< "set cbrange [0:max_velocity]\n\n"
<< "# Set the title for the color palette\n"
<< "set cblabel \"V Magnitude\"\n\n"
<< "# Set the color palette formula\n"
<< "set palette model CMY rgbformulae 7,5,15\n\n"
<< "vector_length = " + std::to_string(vectorLengthMultiplier) + " * max_velocity  # Set your desired constant vector length\n\n"
<< "# Set up the loop for animation\n"
<< "do for [k = 1:" << Nsave+2 << "]{\n"
<< "    plot for [i=0:" << NumberParticles - 1 << "] sprintf('Trajectories/Exercise" + Exercise + "/" + subExercise + "/ST_" + oss.str() + "/data/Trajectory_%d', i) every ::k::k using 2:3 with p pt 7 lc 1 title '', \\\n"
<< "    file_path_flowfield u 1:2:(vector_length * $5 / sqrt($5**2 + $6**2)):(vector_length * $6 / sqrt($5**2 + $6**2)):(sqrt($5**2 + $6**2)) \\\n"
<< "    w vectors filled lc palette title ''\n"
<< "}\n\n"
<< "# Reset terminal and output\n"
<< "set output\n";

// Close the file
plotFile.close();

std::cout << "plotTrajectory.p file created successfully." << std::endl;

// Run gnuplot with the created script
int result2 = system("gnuplot plotTrajectory.p");
//int result2 = system("gnuplot .\\plotTrajectory.p");

// Check if the gnuplot command was successful
if (result2 == 0) {
std::cout << "Gnuplot executed successfully!" << std::endl;
} else {
std::cerr << "Failed to execute gnuplot!" << std::endl;
}
