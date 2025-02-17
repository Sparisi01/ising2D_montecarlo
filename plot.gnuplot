set terminal 'pngcairo'
set output 'file.png'

plot "data.dat" using 1:2 w l t "Magnetizzazione"

plot "data.dat" using 1:3 w l t "Energia"