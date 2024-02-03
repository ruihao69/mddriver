file='lorenz.xyz'
set term gif animate delay 10 size 3200,2400
set output 'lorenz.gif'

set xrange [-20:25]
set yrange [-25:30]
set zrange [0:55]
set xlabel "X"
set ylabel "Y"
set zlabel "Z"

total_frames = 8300  # Total number of frames in the simulation
target_duration = 20  # Target duration of the GIF in seconds
skip_frames = 300  # Number of frames to skip between each plotted frame

# Calculate the number of frames to plot
num_frames = int(total_frames / skip_frames)

# Calculate the frame interval based on the target duration
frame_interval = target_duration * 1000.0 / total_frames

do for [i=0:num_frames-1] {
    frame_index = i * skip_frames

    # plot the entire dynamics in phase space as lines

    # Plot the current frame with points
    splot file u 2:3:4 with lines lw 2.5 lc rgb "blue", \
          file every 30::0::frame_index u 2:3:4 with points pt 7 ps 1.5 lc rgb "black"

    # Pause between frames
    pause frame_interval/1000.0
}

set output
