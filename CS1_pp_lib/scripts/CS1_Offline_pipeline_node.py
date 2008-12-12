from optparse import OptionParser
from LOFAR_Parset import Parset
import sys, os, time

def add_log(options, line):
    logfile = open(options.log, 'a+')
    logfile.write(line + ' ' + time.ctime() + '\n')
    logfile.close()

parser = OptionParser(usage='%prog [options] -m<MS name> -r<remote_host>')
# use when the location of the MS is not mounted locally, or on a slow NFS
parser.add_option("-r", "--remote_host", dest="remote_host", default="",
                  help="Host name")
parser.add_option("-m", "--ms", dest="ms", default="",
                  help="Measurementset name")
parser.add_option("-l", "--log", dest="log", default="",
                  help="Logfile name")
parser.add_option("-d", "--dir", dest="directory", default="",
                  help="Directory name")


options, args = parser.parse_args()

print 'Started processing:'

if not os.path.exists(options.directory):
  os.mkdir(options.directory)
os.chdir(options.directory)

os.system('rm -rf *.MS *.parset *.debug *.log')
add_log(options, 'Cleaned up old data:')

# use when the location of the MS is not mounted locally, or on a slow NFS
#add_log(options, 'scp -r ' + options.remote_host + ':' + options.ms + ' .')
#result = os.system('scp -r ' + options.remote_host + ':' + options.ms + ' .')
#add_log(options, 'Completed copying data:')
#MS=os.path.split(options.ms)[1]
MS=options.ms

add_log(options, "Preparations done")
sys.stdout.flush

print 'Start processing: ' + MS
sys.stdout.flush
add_log(options, 'Start processing: ' + MS)

if MS == "":
    add_log(options, "No Measurement set given, exiting")
    sys.exit(1)
else:
    if not os.path.exists(MS):
        add_log(options, "MS does not exist, exiting")
        sys.exit(2)

try:
    os.environ['AIPSPATH']
except:
    add_log(options, "AIPSPATH not set, exiting")
    sys.exit(3)

#Bandpass_parset = Parset()
#Bandpass_parset['ms']     = MS
#Bandpass_parset['fixed']  = 5
#Bandpass_parset['window'] = 1
#Bandpass_parset.writeToFile("CS1_BandpassCorrector.parset")
#fd = open("CS1_BandpassCorrector.debug", 'w')
#fd.write("Global 20\n")
#fd.write("Everything 20\n")
#fd.close()

#Flagger_parset = Parset()
#Flagger_parset['ms']        = MS
#Flagger_parset['existing']  = False
#Flagger_parset['threshold'] = 2.0
#Flagger_parset.writeToFile("CS1_FrequencyFlagger.parset")
#fd = open("CS1_FrequencyFlagger.debug", 'w')
#fd.write("Global 20\n")
#fd.write("Everything 20\n")
#fd.close()

#Squasher_parset  = Parset()
#(head, tail) = os.path.split(MS)
#Squasher_parset['inms']          = MS
#Squasher_parset['outms']         = tail + "s"
#Squasher_parset['start']         = 0
#Squasher_parset['step']          = 32
#Squasher_parset['nchan']         = 256
#Squasher_parset['threshold']     = 0
#Squasher_parset['useflags']      = True
#Squasher_parset['allcolumns']    = True
#Squasher_parset.writeToFile("CS1_DataSquasher.parset")
#fd = open("CS1_DataSquasher.debug", 'w')
#fd.write("Global 20\n")
#fd.write("Everything 20\n")
#fd.close()

## These are useful imager parameters only around 60 MHz, for other frequencies you need to change some values
#Imager_parset  = Parset()
#Imager_parset['ms']            = MS
#Imager_parset['compress']      = "False"
#Imager_parset['datamode']      = "channel"
#Imager_parset['imagemode']     = "mfs"
#Imager_parset['spwid']         = [0,1]
#Imager_parset['nchan']         = 256
#Imager_parset['start']         = 0
#Imager_parset['step']          = 1
#Imager_parset['nx']            = 512
#Imager_parset['ny']            = 512
#Imager_parset['cellx']         = 750
#Imager_parset['celly']         = 750
#Imager_parset['stokes']        = "I"
#Imager_parset['weighttype']    = "natural"
#Imager_parset['weightnpixels'] = 1024
#Imager_parset['tile']          = 32
#Imager_parset['padding']       = 1.0
#Imager_parset['gridfunction']  = "SF"
#Imager_parset['imagetype']     = "observed"
#Imager_parset['imagename']     = "observed.image"
#Imager_parset.writeToFile("CS1_Imager.parset")
#fd = open("CS1_Imager.debug", 'w')
#fd.write("Global 20\n")
#fd.write("Everything 20\n")
#fd.close()

#os.system("/app/LOFAR/stable/CS1_BandpassCorrector")
#add_log(options, 'CS1_BandpassCorrector finished')
#os.system("/app/LOFAR/stable/CS1_FrequencyFlagger")
#add_log(options, 'CS1_FrequencyFlagger finished')
#os.system("/app/LOFAR/stable/CS1_DataSquasher")
#add_log(options, 'CS1_DataSqasher finished')

IDPPP_parset  = Parset()
(head, tail) = os.path.split(MS)
IDPPP_parset['msin']          = MS
IDPPP_parset['msout']         = tail + "s"
IDPPP_parset['bandpass']      = 1 ##which bandpass to use
IDPPP_parset['flagger']       = 4 ##which flagger to use
IDPPP_parset['squasher']      = 1 ##which squasher to use
IDPPP_parset['fixed']         = 5 ##bandpass selctor
IDPPP_parset['freqwindow']    = 9 ##FrequencyFlagger, MADFlagger
IDPPP_parset['timewindow']    = 7 ##ComplexMedianFlagger(2), MADFlagger
IDPPP_parset['threshold']     = 4 ##FrequencyFlagger, MADFlagger
IDPPP_parset['min']           = 1 ##ComplexMedianFlagger(2)
IDPPP_parset['max']           = 1 ##ComplexMedianFlagger(2)
IDPPP_parset['existing']      = False ##All flaggers
IDPPP_parset['nchan']         = 224 ##Squasher
IDPPP_parset['start']         = 16 ##Squasher
IDPPP_parset['step']          = 64 ##Squasher
IDPPP_parset['skipflags']     = True ##Squasher
IDPPP_parset['allcolumns']    = False ##Squasher
IDPPP_parset.writeToFile("CS1_IDPPP.parset")
fd = open("CS1_IDPPP.debug", 'w')
fd.write("Global 20\n")
fd.write("Everything 20\n")
fd.close()

log_prop = """#
# setup the right levels for logging and tracing
#
log4cplus.rootLogger=DEBUG, STDOUT
log4cplus.logger.TRC=DEBUG, STDOUT
log4cplus.additivity.TRC=FALSE

#
# define the output channels
#
log4cplus.appender.STDOUT=log4cplus::ConsoleAppender
log4cplus.appender.STDOUT.layout=log4cplus::PatternLayout
log4cplus.appender.STDOUT.layout.ConversionPattern=%D{%d-%m-%y %H:%M:%S} %-5p %c{3} - %m [%.25l]%n

#log4cplus.appender.MACSTDERR=log4cplus::ConsoleAppender
#log4cplus.appender.MACSTDERR.layout=log4cplus::PatternLayout
#log4cplus.appender.MACSTDERR.layout.ConversionPattern=%x %D{%d-%m %H:%M:%S.%q} %-5p %c{4} - %m [%.25l]%n
#log4cplus.appender.MACSTDERR.logToStdErr=true

#log4cplus.appender.MACCLP=log4cplus::SocketAppender
#log4cplus.appender.MACCLP.port=23999
#log4cplus.appender.MACCLP.host=mcu001t
#log4cplus.appender.MACCLP.Threshold=INFO

#log4cplus.appender.FILE=log4cplus::RollingFileAppender
#log4cplus.appender.FILE.File=./test.log
#log4cplus.appender.FILE.MaxFileSize=10MB
#log4cplus.appender.FILE.MaxBackupIndex=9
#log4cplus.appender.FILE.layout=log4cplus::PatternLayout
#log4cplus.appender.FILE.layout.ConversionPattern=%x %D{%d-%m-%y %H:%M:%S} %-5p %c{3} - %m [%.25l]%n

#log4cplus.appender.DUMP=log4cplus::NullAppender"""

fd = open("CS1_IDPPP.log_prop", 'w')
fd.write(log_prop)
fd.close()

os.system("/app/LOFAR/stable/CS1_IDPPP")
add_log(options, 'CS1_IDPPP finished')

# Doesn't work on 64 bit systems
#os.system("glish -l /app/LOFAR/stable/flag_auto.g " + MS)
#add_log(options, "Flagging auto correlations done")

#Messes up if we don't flag the auto correlations
#os.system("/app/LOFAR/stable/CS1_Imager")
#add_log(options, "Imager done")

# Not using this when testing
#os.system('rm -rf ' +  MS)
#add_log(options, 'Deleting MS finished')
