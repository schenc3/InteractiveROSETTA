import glob
import os
import os.path
import time
import datetime
import sys
import gzip
from subprocess import Popen, PIPE, STDOUT

# ====================================================================================================================
# These are variables that need to be changed to suit your setup
iRosetta_home = "/bach1/home/schenc3/public_html/InteractiveROSETTA"
# This is not the same file as for the Windows backend
# It is the list of nodes that are MPI-capable for running MPI C++ Rosetta
# mpi_hostlist is simply the hostlist file used for mpirun
mpi_hostlist = iRosetta_home + "/hostbach"
# This value is true if the server running daemon_server.py is not able to submit MPI jobs to a cluster
# If True, you need to specify the name of the machine that can submit MPI jobs so this script can ssh into it
separateMPI_master = True
MPI_master = "bach1"
# MPI-specific commands
# Change if you are not using mpirun or mpirun does not take these arguments
mpiexec = "mpirun --mca btl_tcp_if_exclude lo,virbr0"
hostfile_arg = "--hostfile"
numproc_arg = "-np"
# Location of Rosetta executables
rosetta_bin = "/bach1/usr/local/rosetta_2014.22.56873_bundle/main/source/bin"
# ====================================================================================================================

# This script is used by InteractiveROSETTA to submit uploaded jobs to the C++ Rosetta, and then packages them up into
# files that the client GUI can access remotely
# Please note that this script is called by daemon_server.py
# It is running from iRosetta_home/results/jobID, so iRosetta_home is "../.." from the current directory
if (len(sys.argv) != 3):
    print "Usage: python rosetta_submit.py protocol jobID"
    exit()
# Get the backend hosts
f = open(mpi_hostlist, "r")
hosts = []
for aline in f:
    if (aline.find("#") >= 0):
        continue
    host = aline.split()[0]
    cpus = int(aline.split()[1].split("=")[1])
    hosts.append([host, cpus])
f.close()
protocol = sys.argv[1].strip()
jobID = sys.argv[2].strip()
# ===========================================================================================================================
# How many nodes do we need?  Find out by counting the number of lines in the states files
# Files should have already all been unpacked and we should be running from the jobID result folder already
cpus = 0
if (protocol == "msd"):
    states = glob.glob("*.states")
    for statefile in states:
        f = open(statefile.strip(), "r")
        for aline in f:
            if (len(aline.strip()) > 0):
                cpus = cpus + 1
        f.close()
runnable = False
selected_hosts = []
cpus_allocated = 0
# Update the flags file with the params files
# First check to see if they are already in the flags file (if you have to restart for some reason
# it will add duplicate lines and crash Rosetta)
addparams = True
f = open("flags", "r")
for aline in f:
    if ("-extra_res_fa" in aline):
        addparams = False
        break
f.close()
if (addparams):
    paramsfiles = glob.glob("*.params")
    if (len(paramsfiles) > 0):
        f = open("flags", "a")
        f.write("-extra_res_fa")
        for params in paramsfiles:
            f.write(" " + params.strip())
        f.write("\n")
        f.close()
# Enter the queue
while (os.path.isfile(iRosetta_home + "/writing_to_queue")):
    time.sleep(10)
f = open(iRosetta_home + "/writing_to_queue", "w")
f.write("TAKEN")
f.close()
f = open(iRosetta_home + "/rosetta_queue", "a")
f.write(jobID.strip() + "\n")
f.close()
os.remove(iRosetta_home + "/writing_to_queue")
while (True):
    # We'll use a common file to queue jobIDs so they all run in order
    # First grab a mutex in case multiple instances of this script are running, so they don't query the backends at the same
    # time and grab the same nodes
    while (os.path.isfile(iRosetta_home + "/writing_to_queue")):
        time.sleep(10)
    f = open(iRosetta_home + "/writing_to_queue", "w")
    f.write("TAKEN")
    f.close()
    f = open(iRosetta_home + "/rosetta_queue", "r")
    data = f.readlines()
    if (len(data) > 0 and data[0].strip() == jobID.strip()):
        runnable = True
    else:
        runnable = False
    f.close()
    if (runnable):
        f = open(iRosetta_home + "/rosetta_queue", "w")
        for aline in data[1:]:
            f.write(aline.strip() + "\n")
        f.close()
        os.remove(iRosetta_home + "/writing_to_queue")
        break
    os.remove(iRosetta_home + "/writing_to_queue")
    time.sleep(10)
# Loop for finding free nodes
while (True):
    while (os.path.isfile(iRosetta_home + "/backend_query")):
        time.sleep(10)
    f = open(iRosetta_home + "/backend_query", "w")
    f.write("TAKEN")
    f.close()
    # Query the backend nodes to find which nodes are free and submit to those nodes
    for host in hosts:
        hostname = host[0]
        slots = host[1]
        if (separateMPI_master):
            process = Popen(args="ssh " + MPI_master + " \'ssh " + hostname + " \"uptime\"\'", stdout=PIPE, shell=True)
        else:
            process = Popen(args="ssh " + hostname + " \"uptime\"", stdout=PIPE, shell=True)
        uptimeoutput = process.communicate()[0]
        try:
            load_1min = float(uptimeoutput.split()[len(uptimeoutput.split())-3].split(",")[0])
        except:
            load_1min = 0
        if (float(load_1min) < 0.25):
            cpus_allocated = cpus_allocated + slots
            selected_hosts.append([hostname, slots])
        if (cpus_allocated >= cpus):
            runnable = True
            break
        else:
            runnable = False
    os.remove(iRosetta_home + "/backend_query")
    if (runnable):
        break
    time.sleep(10)
# Generate a hostfile for this job
# This assumes that the hostfile is formatted as: machinename slots=n max_slots=p
# Where machinename is the name of the computer, n is the number of CPUs, and m is the maximum number of processes
# Usually n = m
# If this is not how hostfiles are formatted on your machine, you need to change this code
f = open("hostfile_" + jobID.strip(), "w")
for host in selected_hosts:
    f.write(host[0] + " slots=" + str(host[1]) + " max_slots=" + str(host[1]) + "\n")
f.close()
# End potential hostfile code change region
if (runnable):
    if (protocol == "msd"):
        if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/mpi_msd.mpi.linuxgccrelease @flags > msd.out) >& msd.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/mpi_msd.mpi.linuxgccrelease @flags > msd.out) >& msd.err"
    print command
    p = Popen(args=command, stdout=PIPE, shell=True)
    (out, err) = p.communicate()
# Now package the results up
if (protocol == "msd"):
    # Did we crash?  The answer is yes if msd.err has content in it
    crashed = False
    errmsg = ""
    f = open("msd.err", "r")
    for aline in f:
        if (len(aline.strip()) > 0):
            crashed = True
        errmsg = errmsg + aline.strip() + "\n"
    f.close()
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
        exit()
    outputfiles = glob.glob("*msd_output*.pdb")
    # Rename these files, I hate the standard output names
    # Figure out what the actual PDB name is, and generate the name based on that
    # The Rosetta name is msd_output_{rank}_vec{state}_{indx}.pdb
    # I want {pdbname}_MSD_{rank}.pdb
    for pdbfile in outputfiles:
        rank = pdbfile.split("_")[2]
        vec = pdbfile.split("_")[3]
        indx = int(pdbfile.split("_")[4].split(".pdb")[0])
        # Get the PDB name from the states file
        currindx = 1
        pdbname = "unknown"
        statefile = vec[3:] + ".states"
        f = open(statefile, "r")
        for aline in f:
            if (".pdb" in aline):
                if (currindx == indx):
                    pdbname = aline.split()[0].split(".pdb")[0]
                    break
                else:
                    currindx = currindx + 1
        f.close()
        os.rename(pdbfile, pdbname + "_MSD_" + rank + ".pdb")
    outputfiles = glob.glob("*_MSD_*.pdb")
archive = gzip.open("resultstemp.gz", "wb")
for pdbfile in outputfiles:
    archive.write("BEGIN PDB " + pdbfile.strip() + "\n")
    f = open(pdbfile.strip(), "r")
    for aline in f:
        archive.write(aline.strip() + "\n")
    f.close()
    archive.write("END PDB " + pdbfile.strip() + "\n")
archive.close()
os.rename("resultstemp.gz", "results.gz")
