set term 'cairolatex' pdf dashed lw 2 size 6.8 in, 3.9 in
set output 'survival-plot.tex'

set lmargin 8.1
set rmargin 8

set xtics add (1)
set y2tics out
set logscale y
set format y '#@%.0E'
set format y2 '%4.2f'

set grid

set xlabel '$\sqrt{s}$, GeV'
set ylabel '$\mathrm{d} L / \mathrm{d} \sqrt{s},~\text{GeV}^{-1}$' \
           offset character 3
set y2label 'Survival factor' offset character -1

set key Left reverse samplen 4

plot [1 : 1e3] \
     'luminosity.txt' w l lc rgb 'red' dt (20, 5) \
     t '$\mathrm{d} L / \mathrm{d} \sqrt{s} \rvert_{P = 1}$', \
     'luminosity_b.txt' w l lc rgb 'blue' dt (20, 5, 5, 5) \
     t '$\mathrm{d} L / \mathrm{d} \sqrt{s}$', \
     'parallel.txt' w l lc rgb 'dark-green' dt (20, 20) \
     t '$\mathrm{d} L_\parallel / \mathrm{d} \sqrt{s}$', \
     'perpendicular.txt' w l lc rgb 'brown' dt (5, 5) \
     t '$\mathrm{d} L_\perp / \mathrm{d} \sqrt{s}$', \
     'survival.txt' axes x1y2 w l lc rgb 'black' \
     t 'Survival factor'
