
#          Copyright Balazs Cziraki 2016
# Distributed under the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

reset

set term wxt enh

unset parametric

datafile = "test.txt"

set autoscale xy
set autoscale z
set cbrange[-pi:pi]

load "settings.plt"

set xtics format "%g"
set ytics format "%g"
set ztics format "%g"

set size square
set ticslevel 0

set pm3d corners2color geomean scansautomatic
set palette model HSV defined ( 0 1 1 1 , 1 0 1 1 )
set grid

set sample 121
set iso 121

set xlabel "Re(z)"
set ylabel "Im(z)"
set zlabel "|w(z)|"
set cblabel "Arg(w(z))"

set cbtics ("-{/Symbol p}" -pi, "-{/Symbol p}/2" -pi/2, "0" 0, "{/Symbol p}/2" pi/2, "{/Symbol p}" pi)

i={0,1}
z(x,y)=x+i*y
crd(z)=z
w(z)=z
wi(z)=z
absw(x,y)=abs(w(crd(z(x,y))))
argw(x,y)=arg(w(crd(z(x,y))))
rez(x,y)=real(wi(crd(z(x,y))))
imz(x,y)=imag(wi(crd(z(x,y))))

#load "my_functions.plt"

splot datafile u 1:2:(abs(z($3,$4))):(arg(z($3,$4))) title "w(z)" w pm3d lt 1 lc 0 lw 1

pause -1
