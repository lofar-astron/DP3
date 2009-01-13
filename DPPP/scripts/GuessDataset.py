import sys
print "Script for generating DPPP config file"
sys.stdout.write("Give the number of (micro)stations [16]: ")
stations = sys.stdin.readline()
if len(stations) <= 1: stations = 16
else: stations = int(stations)

sys.stdout.write("Give the integration time [10]: ")
inttime = sys.stdin.readline()
if len(inttime) <= 1: inttime = 10
else: intttime = int(inttime)

sys.stdout.write("Give the number of bands [36]: ")
bands = sys.stdin.readline()
if len(bands) <= 1: bands = 36
else: bands = int(bands)

sys.stdout.write("Give the duration of the measurement [1:23] hh:mm: ")
duration = sys.stdin.readline()
if len(duration) <= 1: duration = "1:23\n"

baselines = stations *(stations+1) /2
print "Total datasize is %s" % (256*9*baselines * bands * int(duration[:-4]) *60 * int(duration[-3:-1]) * inttime/60,)
