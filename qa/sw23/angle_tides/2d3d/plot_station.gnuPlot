


mpl_top    = 0.4 #inch  outer top margin, title goes here
mpl_bot    = 0.9 #inch  outer bottom margin, x label goes here
mpl_left   = 0.5 #inch  outer left margin, y label goes here
mpl_right  = 0.1 #inch  outer right margin, y2 label goes here
mpl_height = 1.5 #inch  height of individual plots
mpl_width  = 2.0 #inch  width of individual plots
mpl_dx     = 0.1 #inch  inter-plot horizontal spacing
mpl_dy     = 0.1 #inch  inter-plot vertical spacing
mpl_ny     = 3   #number of rows
mpl_nx     = 2   #number of columns

# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

# placement functions
#   rows are numbered from bottom to top
bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
#   columns are numbered from left to right
left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)

#set terminal postscript eps enhanced color dl 2.0 size xsize,ysize "Helvetica" 28
#set encoding iso_8859_1
#set tics scale 1.5

#set output 'nxm_plot.eps'

set offsets
set autoscale fix
set size 1,1
set nokey

# define x-axis settings for all subplots
#set xrange [-4:4]
#set xlabel ''
#set format x ''
#set xtics pi
#set mxtics 4

set xtics rotate by 60 right

# start plotting
set multiplot


#-----------------------------------------------
# subplot  1-3
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for third row (top)
set tmargin at screen top(3)
set bmargin at screen bot(3)

set title 'SW 2D'

set ylabel "elevation"
set xlabel ""             # no label here
set format y "%-4.2f"
set format x ""

plot 'station2d.out' using 1:2  with line lt -1 lw 1.5, 'station2d.out' using 1:3 pt 11 ps 0.5 lc rgb "blue"

#-----------------------------------------------
# subplot  2-3
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for third row (top)
set tmargin at screen top(3)
set bmargin at screen bot(3)

set title 'SW 3D'

set ylabel ""             # no label here
set xlabel ""             # no label here
set format y ""           # no tic labels

plot 'station3d.out' using 1:2 with line lt -1 lw 1.5, 'station3d.out' using 1:3 pt 11 ps 0.5 lc rgb "blue"

#-----------------------------------------------
# subplot  1-2

#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)

#  set horizontal margins for second row (middle)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set title ''

set ylabel "u"
set xlabel ""             # no label here

#set yrange [-1.5:1.5]
set format y "%-1.1f"

plot 'station2d.out' using 1:4 with line lt -1 lw 1.5, 'station2d.out' using 1:5 pt 11 ps 0.5 lc rgb "red"

#-----------------------------------------------
# subplot  2-2
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for second row (middle)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set title ''

set ylabel ""             # no label here
set xlabel ""             # no label here
#set yrange [-1.5:1.5]
set format y ""           # no tic labels

plot 'station3d.out' using 1:4 with line lt -1 lw 1.5, 'station3d.out' using 1:5 pt 11 ps 0.5 lc rgb "red"

#-----------------------------------------------

set xlabel

# subplot  1-2
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for first row (bottom)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set title ''

# now set a label and tic marks for the x-axis
set xlabel "time (s)"
set ylabel "v"
set format y "%-1.1f"
set format x "%-4.1e"

plot 'station2d.out' using 1:6 with line lt -1 lw 1.5, 'station2d.out' using 1:7 pt 11 ps 0.5 lc rgb "green"

#-----------------------------------------------
# subplot  2-2
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for first row (bottom)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set title ''

set xlabel "time (s)"
set ylabel ""             # no label here
set format y ""           # no tic labels
set format x "%-4.1e"

plot 'station3d.out' using 1:6 with line lt -1 lw 1.5, 'station3d.out' using 1:7 pt 11 ps 0.5 lc rgb "green"

unset multiplot