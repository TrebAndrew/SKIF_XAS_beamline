# Use this file with GnuPlot to plot the data.
# ATTENTION: GnuPlot v.5 or later is required!
# GnuPlot home is: http://gnuplot.sourceforge.net
set terminal png enhanced small size 448,336 background "#000000"
set output "gd787132.png"
set border 15 lc rgb "#00ff00"
set tics out
set grid lc rgb "#a1a100"
unset key
set style data lines
set format y "%g"
set xlabel "Scan angle [sec]"
set ylabel "Reflectivity"
set logscale y
set yrange [0.0001 :    1]
plot "gd787132.dat" using 1:2 lc rgb "#ff0000"
