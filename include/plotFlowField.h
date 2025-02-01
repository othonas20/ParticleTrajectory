//
// Created by VS0121 on 22/8/2024.
//

// Create and open the cylinderfield.p file
std::ofstream plotFile2(pngName + ".p");

// Check if the file was created successfully
if (!plotFile2) {
std::cerr << "Error creating wallFlow.p file!" << std::endl;
return 1;
}

// Write the provided gnuplot script into the file
plotFile2 << "# Set up the terminal and output file (png)\n"
<< "set terminal png size 1000,800\n"
<< "set output 'Trajectories/Exercise" + Exercise
+ "/" + subExercise
<< + "/ST_" + oss.str() + "/" + pngName + "_Stk_" + oss.str() + ".png'\n\n"
<< "# Set plot properties\n"
<< "set title '" + pngName + " Stk " + oss.str() + "'\n"
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
<< + "/ST_" + oss.str() + "/" + pngName +"_Stk_" + oss.str() + "'\n\n"
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
<< "#set cbrange [0:max_velocity]\n\n"
<< "# Set the title for the color palette\n"
<< "set cblabel \"V Magnitude\"\n\n"
<< "# Set the color palette formula\n"
<< "set palette model CMY rgbformulae 7,5,15\n\n"
<< "vector_length = " + std::to_string(vectorLengthMultiplier) + " * max_velocity  # Set your desired constant vector length\n\n"
<< "plot file_path_flowfield u 1:2:(vector_length * $5 / sqrt($5**2 + $6**2)):(vector_length * $6 / sqrt($5**2 + $6**2)):(sqrt($5**2 + $6**2)) \\\n"
<< "    w vectors filled lc palette title ''\n\n"
<< "# Reset terminal and output\n"
<< "set output\n";

// Close the file
plotFile2.close();

std::cout << pngName + ".p file created successfully." << std::endl;

// Run gnuplot with the created script
//int result3 = system("gnuplot cylinderFlow.p");
int result3 = system("gnuplot vortexFlow.p");
//int result3 = system("gnuplot .\\cylinderFlow.p");

// Check if the gnuplot command was successful
if (result3 == 0) {
std::cout << "Gnuplot executed successfully!" << std::endl;
} else {
std::cerr << "Failed to execute gnuplot!" << std::endl;
}