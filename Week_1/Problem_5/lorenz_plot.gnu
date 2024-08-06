set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'lorenz_plot.png'

set title "Lorenz Attractor"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"

set grid

splot 'lorenz_sol1.txt' using 2:3:4 with lines title 'Initial Condition 1', \
      'lorenz_sol2.txt' using 2:3:4 with lines title 'Initial Condition 2'
