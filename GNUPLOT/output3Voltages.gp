
set terminal pngcairo size 1000,600 enhanced font 'Verdana,10'
set output 'neuron_voltages.png'

set title "Neuron Membrane Voltages Over Time"
set xlabel "Time (ms)"
set ylabel "Voltage (mV)"

set grid
set key outside right top

set datafile separator whitespace

plot 'NeuronVoltages.txt' using 1:2 with lines lw 2 title 'Neuron 1', \
     ''         using 1:3 with lines lw 2 title 'Neuron 2', \
     ''         using 1:4 with lines lw 2 title 'Neuron 3'