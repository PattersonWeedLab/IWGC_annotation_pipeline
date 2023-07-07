#!/bin/bash
LOG_FILE="$2"
PID="$1"

top -b -d 2 -p $PID | awk \
    -v cpuLog="$LOG_FILE" -v pid="$PID" -v pname="$PNAME" '
    /^top -/{time = $3}
    $1+0>0 {printf "%s %s :: %s[%s] CPU Usage: %d%%\n", \
            strftime("%Y-%m-%d"), time, pname, pid, $9 > cpuLog
            fflush(cpuLog)}'
