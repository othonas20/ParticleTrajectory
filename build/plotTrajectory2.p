# Set up the terminal and output file (GIF)
set terminal gif animate optimize delay 10 size 2000,2000
set output 'Trajectories/ExerciseC/3/ST_10/trajectory_animation_C_ST_10.gif'

# Set plot properties
set title "Trajectory Stk = 10"
set xlabel "X"
set ylabel "Y"
set grid

set size ratio -1
unset colorbox  # Remove the default colorbox (used for vectors)

set pointsize 1
set style line 2 lc rgb '#0060ad' pt 7

# Define the animation parameters
file_path_flowfield = 'Trajectories/ExerciseC/3/ST_10/cylinderFlow_Stk_10'

# Calculate the maximum magnitudes for vectors
stats file_path_flowfield using (sqrt($5**2 + $6**2)) name "Velocities"
stats file_path_flowfield using 1 name "X"
stats file_path_flowfield using 2 name "Y"

max_velocity = Velocities_max
DX = X_max - X_min
DY = Y_max - Y_min

set xrange [X_min - 0.2*DX : X_max + 0.2 * DX]
set yrange [Y_min - 0.2 * DY : Y_max + 0.2 * DY]

# Set the color palette for the temperature field
set cbrange [300:900*0.7]  # Adjust temperature range as needed
set cblabel "Temperature (K)"
set colorbox

# Set color palette for temperature (surface)
#set palette rgbformulae 34,35,36
set palette defined (0 'white', 0.5 'yellow', 1 'orange')

# Set a fixed vector length
vector_length = 0.000050 * max_velocity

# Set up the loop for animation
do for [k = 1:202] {
    # Plot the surface temperature field using column 1:2:8 from the data file
    plot file_path_flowfield u 1:2:8 w image notitle, \
    for [i=0:49] sprintf('Trajectories/ExerciseC/3/ST_10/data/Trajectory_%d', i) every ::k::k u 2:3 w p pt 7 lc 1 title '', \
    file_path_flowfield u 1:2:(vector_length * $5 / sqrt($5**2 + $6**2)):(vector_length * $6 / sqrt($5**2 + $6**2)):(sqrt($5**2 + $6**2)) w vectors lc 'black' nohead title ''
}

# Reset terminal and output
set output
