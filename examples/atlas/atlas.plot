set term 'cairolatex' dashed lw 2 size 6.8 in, 4.0 in
set output 'atlas-plot.tex'

set multiplot

s = 0.3
set size 1, 1-s
set origin 0, s

set lmargin 7
set bmargin 0

set format x ''
set format y '%.2f'

set ylabel '$\mathrm{d} \sigma_\text{fid.}(pp \to pp \mu^+ \mu^-) / \mathrm{d} \sqrt{s}$, pb/GeV' offset 0.9, 1.2

set grid

set key Left reverse

plot [10:] \
     'experiment.txt'     u 1:($2*1e12):($3*1e12) w yerror \
     t 'Experimental data', \
     'differential.txt'   u 1:($2*1e12) w l lc rgb 'red' dt (20, 5) \
     t 'Non-electromagnetic interactions are neglected', \
     'differential_b.txt' u 1:($2*1e12) w l lc rgb 'black' \
     t 'Non-electromagnetic interactions are respected', \
     'integrated.txt'     u 1:($2*1e12) w steps lc rgb 'red' dt (20, 5) notitle, \
     'integrated_b.txt'   u 1:($2*1e12) w steps lc rgb 'black' notitle

set size 1, s
set origin 0, 0

set tmargin 0
set bmargin 3

set xlabel '$\sqrt{s}$, GeV'
set ylabel 'EPA / ATLAS' offset 0

set xtics 20, 10
set xtics add (12)
set ytics 0.8, 0.1, 1.1

set format x
set format y '%3.1f'

plot [10:] [0.8:1.3] \
     'experiment.txt' u 1:(1):($3/$2) w yerror notitle, \
     'ratio.txt'      u 1:2 w steps lc rgb 'red' dt (20, 5) notitle, \
     'ratio_b.txt'    u 1:2 w steps lc rgb 'black' notitle
