reset
unset border
set polar
set angles degrees #set gnuplot on degrees instead of radians

set style line 10 lt 1 lc 0 lw 0.2 #redefine a new line style for the grid

set grid polar 30 #set the grid to be displayed every 30 degrees
set grid ls 10

set xrange [-16:16] #make gnuplot to go until 14
set yrange [-16:16]

#set xtics axis #display the xtics on the axis instead of on the border
set ytics axis

set xtics scale 0 #"remove" the tics so that only the y tics are displayed
set xtics ("" 0, "" 2, "" 4, "" 6, "" 8, "" 10, "" 12, "" 14) #set the xtics only go from 0 to 6000 with increment of1000 but do not display anything. This has to be done otherwise the grid will not be displayed correctly.
set ytics 0, 2, 14 #make the ytics go from the center (0) to 10 with increment of 2

set size square 

unset key

set_label(x, text) = sprintf("set label '%s' at (16*cos(%f)), (16*sin(%f))     center", text, x, x) #this places a label on the outside

#here all labels are created
eval set_label(0, "0")
eval set_label(30, "30")
eval set_label(60, "60")
eval set_label(90, "90")
eval set_label(120, "120")
eval set_label(150, "150")
eval set_label(180, "180")
eval set_label(210, "210")
eval set_label(240, "240")
eval set_label(270, "270")
eval set_label(300, "300")
eval set_label(330, "330")

set style line 11 lt 1 lw 2 pt 7 ps 1 #set the line style for the plot

#and finally the plot
set term pngcairo
set output 'windfreq.png'
plot 'horns_rev_windrose.dat' u 1:3 w lp ls 11
