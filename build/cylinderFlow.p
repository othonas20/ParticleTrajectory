# Set up the terminal and output file (png)
set terminal png size 2000,2000
set output 'Trajectories/ExerciseC/3/ST_10/cylinderFlow_Stk_10.png'

# Set plot properties
set title 'cylinderFlow Stk 10'
set xlabel "X"
set ylabel "Y"
set grid

set size ratio -1
set colorbox

set pointsize 1
set style line 2 lc rgb '#0060ad' pt 7

# Define the animation parameters
file_path_flowfield = 'Trajectories/ExerciseC/3/ST_10/cylinderFlow_Stk_10'

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
#set cbrange [0:max_velocity]

# Set the title for the color palette
set cblabel "V Magnitude"

# Set the color palette formula
set palette model CMY rgbformulae 7,5,15

vector_length = 0.000050 * max_velocity  # Set your desired constant vector length

plot file_path_flowfield u 1:2:(vector_length * $5 / sqrt($5**2 + $6**2)):(vector_length * $6 / sqrt($5**2 + $6**2)):(sqrt($5**2 + $6**2)) \
    w vectors filled lc palette title ''

# Reset terminal and output
set output
