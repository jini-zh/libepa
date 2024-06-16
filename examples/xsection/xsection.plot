set term 'cairolatex' pdf dashed lw 2 size 6.8 in, 3.7 in
set output 'xsection-plot.tex'

set rmargin 8

set format y  '%3.1f'
set format y2 '%4.2f'

set xlabel '$m_\chi$, GeV'
set ylabel '$\sigma(pp \to pp \chi^+ \chi^-)$, fb'
set y2label 'Ratio' offset -1

set xtics add (90, 250)
set y2tics out

set grid

set key Left reverse

plot [90:250] \
     'xsection.txt' u 1:($2 * 1e15) w l lc rgb 'red' dt (20, 6) \
     t 'Cross section with non-electromagnetic interactions neglected', \
     'xsection_b.txt' u 1:($2 * 1e15) w l lc rgb 'blue' dt (20, 5, 5, 5) \
     t 'Cross section with non-electromagnetic interactions respected', \
     'ratio.txt' w l axes x1y2 lc rgb 'black' \
     t 'Ratio'
