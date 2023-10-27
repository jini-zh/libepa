set term 'cairolatex' dashed lw 2 size 6.8 in, 3.7 in
set output 'integration-plot.tex'

set format y  '%3.1f'

set xlabel '$a$'
set ylabel '$I(a)$'

set grid

set key Left reverse

set samples 1000

plot [0:100] [0:1.2] \
     x * (x+0.5) / (x+1)**2 lc rgb 'black' t 'Analytical', \
     'integration.txt' pt 1 ps 0.4 lc rgb 'blue' t 'Numerical'
