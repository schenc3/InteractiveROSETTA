#!/bin/csh

# This script should be run as a cronjob every minute or so

# First we have to check to make sure the daemons are running on each of the hosts
# The list of hosts is in hostlist (if you want to add or remove hosts)
# The daemons are programmed to terminate after running for a day so their memory doesn't get stale
# That can be changed if we want

# ==============================================================================================================
# These variables need to be set to their proper values
# This script expects HOSTFILE to be a file containing a list of the names of all the machines that can act as Windows backend
set HOSTFILE = "hostlist"
set REMOTE_SERVER = 0 # Set this to zero if the machine executing this script is the server, 1 if this machine will tell another machine to run the server                                                                                                         
set IROSETTA_HOME = `pwd` # Location of InteractiveROSETTA on this machine                                                        
set REMOTE_IROSETTA_HOME = `pwd` # This is the location on the remote machine, if you are cross mounting from a main file server this needs to be set properly
# ==============================================================================================================
foreach HOST (`cat $HOSTFILE`)
    if ($REMOTE_SERVER == 0) then
        set DAEMON_ON = `ps aux | grep daemon_server | grep -v grep | wc -l` # Should be 1 if the daemon is there                 
        if ($DAEMON_ON == 0) then
            # Not running, restart it                                                                                             
            # Output is redirected to /dev/null because if it is saved to a log file the log file will get big fast               
            # If you need to look at the output because of errors, turn this file off by putting exit at the top of the script    
            # and then run the daemon as a user without redirecting stdout to anything so you can watch the output for specific jobs                                                                                                                               
            cd $IROSETTA_HOME; python daemon_server.py >& $IROSETTA_HOME/daemon_log &
        endif
    else
        set DAEMON_ON = `ssh $HOST "ps aux | grep daemon_server | grep -v grep | wc -l"` # Should be 1 if the daemon is there     
        if ($DAEMON_ON == 0) then
            # Not running, restart it                                                                                             
            # Output is redirected to /dev/null because if it is saved to a log file the log file will get big fast               
            # If you need to look at the output because of errors, turn this file off by putting exit at the top of the script    
            # and then run the daemon as a user without redirecting stdout to anything so you can watch the output for specific jobs                                                                                                                               
            ssh $HOST "cd $REMOTE_IROSETTA_HOME; python daemon_server.py >& $IROSETTA_HOME/daemon_log &"
        endif
    endif
end

# Now we have to clean up all those output files that get generated, because they can build up fast
# If they are more than 5 minutes old, then remove them
# The GUI clients should be able to pick them up within 5 minutes
# If they don't, then the user must have lost their Internet connection in which case the GUI should have decided to use
# a local daemon so they'll get their results locally
foreach OLDFILE (`find $IROSETTA_HOME/results/* -maxdepth 0 -mmin +5 -type f`)
    rm -f $OLDFILE
end
# Sometimes the mutexes get abandoned, and then the whole server hangs until the sysadmin 
# manually deletes them.  This deletes mutexes more than 10 minutes old
foreach OLDFILE (`find $IROSETTA_HOME/backend_query -maxdepth 0 -mmin +10 -type f`)
    rm -f $OLDFILE
end
foreach OLDFILE (`find $IROSETTA_HOME/backend_query -maxdepth 0 -mmin +10 -type f`)
    rm -f $OLDFILE
end
# Clean up some more junk
foreach OLDFILE (`find $IROSETTA_HOME/*.pdb -maxdepth 0 -mmin +20 -type f`)
    rm -f $OLDFILE
end
foreach OLDFILE (`find $IROSETTA_HOME/*-weights -maxdepth 0 -mmin +20 -type f`)
    rm -f $OLDFILE
end
foreach OLDDIR (`find $IROSETTA_HOME/results/*/results.ensb -mtime +7 -type f | awk -F "/results.ensb" '{print $1}'`)
    if (`echo $OLDDIR | grep -c "results"` == 1) then
	rm -fr $OLDDIR
    endif
end
foreach OLDDIR (`find $IROSETTA_HOME/results/*/results.msdar -mtime +7 -type f | awk -F "/results.msdar" '{print $1}'`)
    if (`echo $OLDDIR | grep -c "results"` == 1) then
	rm -fr $OLDDIR
    endif
end
foreach OLDDIR (`find $IROSETTA_HOME/results/*/results.scan -mtime +7 -type f | awk -F "/results.scan" '{print $1}'`)
    if (`echo $OLDDIR | grep -c "results"` == 1) then
	rm -fr $OLDDIR
    endif
end
# Now let's update the statuses of all the jobs that are running
foreach DIR (`ls -d -- results/*/`)
    if (-e $DIR/msd.out) then
	if (`tac $DIR/msd.out | grep -c "Generation"` == 0) then
	    echo "Started MSD" > $DIR/status
	else
	    set GEN = `tac $DIR/msd.out | grep "Generation" | grep -v "Final" | head -n 1 | awk -F "Generation" '{print $2}' | awk '{print $1}'`
	    echo "Generation "$GEN" Completed" > $DIR/status
	endif
    endif
    if (-e $DIR/antibodyinput) then
	if (-e $DIR/antibody.out) then
	    set MODELS = `ls $DIR/*MYAB*.pdb | wc -l`
	    echo $MODELS" Models Completed" > $DIR/status
	else
	    echo "Pending"
	endif
    endif
    if (-e $DIR/coarsedockinput) then
	if (-e $DIR/dock.out) then
	    set R2MODELS = `ls $DIR/*_dock_????_????.pdb | wc -l`
	    if (-e $DIR/round2) then
		echo "Completed" > $DIR/status
	    else if ($R2MODELS == 0) then
		set MODELS = `ls $DIR/*_dock_????.pdb | wc -l`
		echo "Round 1, "$MODELS" Models Completed" > $DIR/status
	    else
		echo "Round 2, "$R2MODELS" Models Completed" > $DIR/status
	    endif
	else
	    echo "Pending"
	endif
    endif
    if (-e $DIR/backrubinput) then
	if (-e $DIR/backrub.out) then
	    set MODELS = `ls $DIR/*.pdb | wc -l`
	    set MODELS = `expr $MODELS - 1`
	    echo $MODELS" Models Completed" > $DIR/status
	else
	    echo "Pending"
	endif
    endif
    if (-e $DIR/flexpepinput) then
	if (-e $DIR/frag.out) then
	    echo "Generating fragments..." > $DIR/status
	else
	    echo "Pending"
	endif
	if (-e $DIR/flexpep.out) then
	    set MODELS = `ls $DIR/*.pdb | wc -l`
	    set MODELS = `expr $MODELS - 1`
	    echo $MODELS" Models Completed" > $DIR/status
	endif
    endif
end
