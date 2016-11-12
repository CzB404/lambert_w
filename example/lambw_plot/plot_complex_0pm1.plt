
#  Copyright Balazs Cziraki 2016.
#  Use, modification and distribution are subject to the
#  Boost Software License, Version 1.0. (See accompanying file
#  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

reset

set term wxt enh
unset surf
set pm3d
set pm3d corners2color geomean
set pm3d depthorder
set palette model HSV defined (0 1 1 1, 1 0 1 1)
set ticslevel 0
set grid

set xlabel "Re(z)"
set ylabel "Im(z)"
set zlabel "|w(z)|"
set cblabel "Arg(w(z))"

set cbrange[-pi:pi]
set cbtics format "%1.1P{/Symbol p}" auto 0.5*pi

i = {0,1}

d=1.0

set xrange[-d:d]
set yrange[-d:d]
set zrange[0:2*d]

set encoding iso_8859_1

set term pngcairo enh truecolor size 1920,1080
set output "lw_composite.png"

set view 70, 130

splot "lw_0.txt" u 1:(-$2):(abs($3+$4*i)):(-arg($3+$4*i)) w pm3d title "W_0(z)",\
"lw_pm1.txt" u 1:2:(abs($3+$4*i)):(arg($3+$4*i)) w pm3d title "W_{\261 1}(z)"

set output
set term wxt enh