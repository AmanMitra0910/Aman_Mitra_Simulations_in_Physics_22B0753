set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'lorenz_plot2.png'

set title "Lorenz Attractor"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"

set grid

splot 'lorenz_sol2.txt' using 2:3:4 with lines title 'Initial Condition 2'
