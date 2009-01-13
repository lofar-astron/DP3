from optparse import OptionParser
from LOFAR_Parset import Parset
import sys, os, time

def add_log(options, line):
    logfile = open(options.run_id + '.log', 'a+')
    logfile.write(line + ' ' + time.ctime() + '\n')
    logfile.close()

def split_remote_path(path):
    # strip trailing white space
    path = path.rstrip()
    # remove trailing slash
    if path[-1] == '/':
        path = path[:-1]
    tmp = path.split(':')
    if len(tmp) == 2:
        return (tmp[0], tmp[1])
    else:
        return ("", path)

def read_file_list(filename):
    flist = []
    dlist = []
    fid = open(filename)
    for line in fid:
        line.strip()
        # skip empty lines
        if len(line) == 0 or line.isspace():
            continue
        if line[0]=='#':
            continue
        tmp = line.split()
        assert len(tmp) == 2
        dlist.append(split_remote_path(tmp[0]))
        flist.append(split_remote_path(tmp[1]))
    fid.close()
    print flist
    print dlist
    return (flist, dlist)

def spawn(dest, src, log):
    try:
        pid = os.fork()
    except OSError, e:
        print 'Unable to fork:' + str(e)
        os._exit(1)
    if pid == 0: ## the new client
        ##Copy node script so we are sure we are running the right version
        os.system('ssh -A -t ' + dest[0] + ' "mkdir '+ dest[1] + '"')
        os.system('scp CS1_Offline_pipeline_node.py ' + dest[0] + ':' + dest[1] +
        '/CS1_Offline_pipeline_node.py')
        ##Run the node script on every node
        if src[0]:
            os.system('ssh -A -t ' + dest[0] + ' "source /app/scripts/doPython;' +
            'source /app/scripts/doStableAIPS++; source /app/scripts/doPyrap; cd '+ dest[1] +
            ';python CS1_Offline_pipeline_node.py -r' + src[0] + ' -m' + src[1] + ' -d' + dest[0] +
            ' -l' + log + ' >> pipeline.log"')
        else:
            os.system('ssh -A -t ' + dest[0] + ' "setenv PYTHONPATH /app/LOFAR/stable;' +
            'source /app/scripts/doStableAIPS++; source /app/scripts/doPyrap; cd '+ dest[1] +
            ';python CS1_Offline_pipeline_node.py  -m' + src[1] + ' -d' + dest[0] +
            ' -l' + log + ' >> pipeline.log"')
## bash version for Pandey
##        os.system('ssh -A -t ' + dest[0] + ' "export PYTHONPATH=/app/LOFAR/stable;' +
##        '. /app/aips++/Unstable/aipsinit.sh; cd '+ dest[1] + ';python ' +
##        'CS1_Offline_pipeline_node.py -r ' + src[0] + ' -m ' + src[1] + ' > pipeline.log"')
        zzz = time.gmtime()
        os.system('scp ' + dest[0] + ':' + dest[1] + '/' + dest[0] +
        '/' + log + ' ' + dest[0] +'_' + str(zzz[3]) + '-' + str(zzz[4]) +
        '-' + str(zzz[5]) + '.log')
        os._exit(0)
    else: ## parent process
        return pid

parser = OptionParser(usage='%prog [options] -m<Obs name>')
parser.add_option("-f", "--file_list", dest="file_list", default="",
                  help="Name of file list")
parser.add_option("-r", "--run_id", dest="run_id", default="",
                  help="ID to identify run")
options, args = parser.parse_args()

try:
    os.environ['AIPSPATH']
except:
    print "AIPSPATH not set, exiting"
    sys.exit(3)


if not ('SSH_AGENT_PID' in os.environ or 'SSH_AUTH_SOCK' in os.environ):
    print "SSH-AGENT not active, exiting"
    sys.exit(3)

run_id = options.run_id

add_log(options, 'Started processing')

(flist, dlist) = read_file_list(options.file_list)

joblist  = []

for i in range(0, len(flist)):
    joblist.append((dlist[i][0], dlist[i][1], flist[i][0], flist[i][1]))


while len(joblist):
    pidlist  = []
    proclist = []
    for job in joblist:
        if not job[0] in proclist:
            pidlist.append((spawn((job[0], job[1]),(job[2],job[3]), run_id + '.log'), job))
            proclist.append(job[0])
        else:
            print "queueing for later execution: " + str(job)
        time.sleep(1) ## If we start all node scripts at the same time, we run into DoS mitigation measures stopping us.
    #
    time.sleep(5) ## wait for childs to spawn, maybe a niceer way to do this?
    #
    for p in pidlist: ## wait for the clients to finish
        (pid, status) = os.waitpid(p[0], 0)
        joblist.remove(p[1])
        print 'result: ' + str(pid) + ' ' + str(status)

add_log(options, "Done with the collapsing of individual subbands")
print "Done with the collapsing of individual subbands "
sys.stdout.flush

add_log(options, "Now copying back The collapsed individual subbands ")
print "Now copying back The collapsed individual subbands "
sys.stdout.flush


MSlist = []
for i in range(0,len(flist)):
    ## get the squashed measurementsets back
    ret = os.system('scp -r ' + dlist[i][0] +':'  +  dlist[i][1] + '/' +
                    dlist[i][0] + '/' + flist[i][1].split('/')[-1]  + 's .')
    print ret
    MSlist.append(flist[i][1].split('/')[-1] + 's')

MS = 'result.MS'

add_log(options, "The copying back is over, preparing to run last steps")
print "The copying back is over, preparing to run last steps "
sys.stdout.flush

Combiner_parset = Parset()
Combiner_parset['inms']  = MSlist
Combiner_parset['outms']   = MS
Combiner_parset.writeToFile("CS1_SPWCombine.parset")
fd = open("CS1_SPWCombine.debug", 'w')
fd.write("Global 20\n")
fd.write("Everything 20\n")
fd.close()

# The old BBS is no longer available. Needs replacing.
#BBS_parset  = Parset()
#BBS_parset['user']                = 'renting'
#BBS_parset['instrumentPDB']       = '/dop58_0/3C343.instr_model'
#BBS_parset['skyPDB']              = '/dop58_0/3C343.sky_model'
#BBS_parset['measurementSet']      = '/dop58_0/CS1-3_12h_256fr.MS'
#BBS_parset['instrumentModelType'] = 'DIPOLE' ##Other options are TOTALEJ,PATCHEJ,REALIMAG,DIPOLE,BANDPASS, sytax for multiple is: DIPOLE.BANDPASS, also TOTALGAIN seems to exist
#BBS_parset['calcUVW']             = false ##this crashes BBS if set to true
#BBS_parset['operation']           = 'solve'
#BBS_parset['columnNameIn']        = 'CORRECTED_DATA'
#BBS_parset['columnNameOut']       = 'MODEL_DATA'
#BBS_parset['timeDomainSize']      = 60
#BBS_parset['startChan']           = 0
#BBS_parset['endChan']             = 255
#BBS_parset['solvParms']           = '[gain:*]'
#BBS_parset['exclParms']           = '[phase:RT0]'
#BBS_parset['nrSolveInterval']     = '[1,6]' ##frequency domain subdivision, time domain subdivision
#BBS_parset['nriter']              = 10
#BBS_parset['saveSolution']        = true


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
#
add_log(options, "Running the combiner")
print "Running the combiner "
sys.stdout.flush

## some processing on the final MS
os.system("/app/LOFAR/stable/CS1_SPWCombine")
add_log(options, "Combiner done")
print "Combiner done "
sys.stdout.flush

#Not working on 64 bit systems
#os.system("glish -l /app/LOFAR/stable/flag_auto.g " + MS)
#add_log(options, "Flagging auto correlations done")
#print "Flagging auto correlations done "
#sys.stdout.flush

#Messes up if we don't flag the auto correlations
#os.system("/app/LOFAR/stable/CS1_Imager")
#add_log("Imager done")
#print "Imager done "
#sys.stdout.flush
