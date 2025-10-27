set terminal png size 800,600
set output 'FrequencyOutput.png'

set xlabel "gSR mS/cm^2"
set ylabel "Frequency (Hz)"

set xrange [.20:.30]
set yrange  [0:9] 

plot "frequency.txt" w l lt rgb "#0000"