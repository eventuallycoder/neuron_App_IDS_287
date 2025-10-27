set terminal png size 800,600
set output 'voltageOutput.png'

set xlabel "time (ms)"
set ylabel "voltage (mV)"

set xrange [0:20000] 

plot "voltage.txt" w l lt rgb "#0000"