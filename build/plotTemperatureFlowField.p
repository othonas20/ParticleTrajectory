set terminal pngcairo size 800,600 enhanced font 'Arial, 12'
set output 'Trajectories/ExerciseC/2/ST_0.5/temperature_field.png'
set size ratio -1

# Set labels for axes
set xlabel 'x'
set ylabel 'y'

# Set color palette for the temperature field
set palette defined (300 "blue", 310 "green", 320 "yellow", 330 "orange", 340 "red")

# Use a color bar to represent the temperature field
set colorbox

# Plot the data as a heatmap for temperature (T)
plot 'Trajectories/ExerciseC/2/ST_0.5/cylinderFlow_Stk_0.5' using 1:2:8 with image notitle
