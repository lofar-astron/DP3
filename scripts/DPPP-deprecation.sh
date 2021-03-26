#!/bin/bash
echo ===================================================================================
echo Warning: you are using DPPP, the deprecated name for DP3, to start the program!
echo In March 2021, the main DP3 executable was renamed to \'DP3\' to use one consistent
echo name. Your parameters are now forwarded to DP3, but this forwarding wrapper will be
echo removed in the future: update your commands/scripts and change \'DPPP\' to \'DP3\'!
echo ===================================================================================
scriptpath=`dirname "$0"`
${scriptpath}/DP3 "$@"
