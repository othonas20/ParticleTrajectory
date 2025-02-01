# Set up the terminal and output file (GIF)
set terminal gif animate optimize delay 10 size 2000,2000
set output 'Trajectories/ExerciseA/1/ST_0.0/trajectory_animation_C_ST_0.0.gif'

# Set plot properties
set title "Trajectory Stk = 0.0"
set xlabel "X"
set ylabel "Y"
set grid

set size ratio -1
set colorbox

set pointsize 1
set style line 2 lc rgb '#0060ad' pt 7

# Define the animation parameters
file_path_flowfield = 'Trajectories/ExerciseA/1/ST_0.0/vortexFlow_Stk_0.0'

# Calculate the maximum magnitudes
stats file_path_flowfield using (sqrt($5**2 + $6**2)) name "Velocities"
stats file_path_flowfield using 1 name "X"
stats file_path_flowfield using 2 name "Y"

max_velocity = Velocities_max
DX = X_max - X_min
DY = Y_max - Y_min

set xrange [X_min - 0.2*DX : X_max + 0.2 * DX]
set yrange [Y_min - 0.2 * DY : Y_max + 0.2 * DY]

# Set the color palette range based on the maximum magnitude
set cbrange [0:max_velocity]

# Set the title for the color palette
set cblabel "V Magnitude"

# Set the color palette formula
set palette model CMY rgbformulae 7,5,15

vector_length = 5.000000 * max_velocity  # Set your desired constant vector length

# Set up the loop for animation
do for [k = 1:202]{
    plot for [i=0:99] sprintf('Trajectories/ExerciseA/1/ST_0.0/data/Trajectory_%d', i) every ::k::k using 2:3 with p pt 7 lc 1 title '', \
    file_path_flowfield u 1:2:(vector_length * $5 / sqrt($5**2 + $6**2)):(vector_length * $6 / sqrt($5**2 + $6**2)):(sqrt($5**2 + $6**2)) \
    w vectors filled lc palette title ''
}

# Reset terminal and output
set output
