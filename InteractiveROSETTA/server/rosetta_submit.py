import glob
import os
import os.path
import shutil
import time
import datetime
import sys
import gzip
import socket
import commands
from subprocess import Popen, PIPE, STDOUT

# ====================================================================================================================
# These are variables that need to be changed to suit your setup
iRosetta_home = "/home/balto/server"
# This is not the same file as for the Windows backend
# It is the list of nodes that are MPI-capable for running MPI C++ Rosetta
# mpi_hostlist is simply the hostlist file used for mpirun
mpi_hostlist = iRosetta_home + "/hostfile"
# This value is true if the server running daemon_server.py is not able to submit MPI jobs to a cluster
# If True, you need to specify the name of the machine that can submit MPI jobs so this script can ssh into it
separateMPI_master = False
MPI_master = "bach1"
this_host = socket.gethostname()
# MPI-specific commands
# Change if you are not using mpirun or mpirun does not take these arguments
mpiexec = "mpirun --mca btl_tcp_if_exclude lo,virbr0"
hostfile_arg = "--hostfile"
numproc_arg = "-np"
# The minimum number of CPUs for parallel jobs that distribute jobs to many CPUs
mincpus = 5
# Location of Rosetta executables
rosetta_bin = "/usr/local/rosetta/main/source/bin"
rosetta_db = "/usr/local/rosetta/main/database"
# Location of the antibody.py script and the antibody/blast databases
antibody_home = "/usr/local/rosetta/tools/antibody"
# ====================================================================================================================

def AA3to1(resn):
    indx3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR HOH ".find(resn)
    if (indx3 < 0):
	return "Z"
    else:
	indx = indx3 / 4
	return "ACDEFGHIKLMNPQRSTVWYO"[indx]

def twoStageDock(nproc, inputstruct, hostfile, receptorchains, ligandchains, decoys, nfinal, randomize1, randomize2, ensemble1, ensemble2):
    outputdir = "."
    # If either ensemble1 or ensemble2 is present, then we need a file for both ensembles
    # If either is False, then that ensemble needs a file with only the single template from the input 
    # structure in it
    if (ensemble1 and not(ensemble2)):
	fout = open("ligand.pdb", "w")
	fin = open(inputstruct, "r")
	for aline in fin:
	    if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] in ligandchains):
		fout.write(aline)
	fin.close()
	fout.close()
	fout = open("ensemble2", "w")
	fout.write("ligand.pdb\n")
	fout.close()
    elif (ensemble2 and not(ensemble1)):
	fout = open("receptor.pdb", "w")
	fin = open(inputstruct, "r")
	for aline in fin:
	    if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] in receptorchains):
		fout.write(aline)
	fin.close()
	fout.close()
	fout = open("ensemble1", "w")
	fout.write("receptor.pdb\n")
	fout.close()
    # Now we have to take some error-preventative measures for ensembles...
    # Since the docked PDB could have come from multiple sources (i.e. the receptor and ligand were in 
    # different PDB files), some chain renaming may have occurred if two chains had the same ID from these
    # files.  Which means the chainIDs don't match what is in the ensemble files
    # So we have to detect this and rename these chains in the ensembles
    if (ensemble1 or ensemble2):
	for i in range(1, 3):
	    if (i == 1):
		inputchains = receptorchains
	    else:
		inputchains = ligandchains
	    # Get the PDBs in the ensemble
	    fin = open("ensemble" + str(i), "r")
	    pdblist = []
	    for aline in fin:
		if (len(aline.strip()) > 0):
		    pdblist.append(aline.strip())
	    fin.close()
	    # Read the chains+sequences in both the input structure and the ensembles
	    # Then we'll assign chains by looking for matches between sequences
	    modeldata = {}
	    ensbdata = {}
	    fin = open(inputstruct, "r")
	    last_res = "00000"
	    for aline in fin:
		if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] in inputchains and aline[22:27] != last_res):
		    chain = aline[21]
		    last_res = aline[22:27]
		    if (len(chain.strip()) == 0):
			chain = "_"
		    try:
			modeldata[chain] += AA3to1(aline[17:20])
		    except:
			modeldata[chain] = AA3to1(aline[17:20])
	    fin.close()
	    fin = open(pdblist[0], "r")
	    last_res = "00000"
	    for aline in fin:
		if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] in inputchains and aline[22:27] != last_res):
		    chain = aline[21]
		    last_res = aline[22:27]
		    if (len(chain.strip()) == 0):
			chain = "_"
		    try:
			ensbdata[chain] += AA3to1(aline[17:20])
		    except:
			ensbdata[chain] = AA3to1(aline[17:20])
	    fin.close()
	    # Now let's figure out the chain mapping
	    chainmap = {}
	    for modelchain in modeldata.keys():
		for ensbchain in ensbdata.keys():
		    if (modeldata[modelchain] == ensbdata[ensbchain]):
			chainmap[ensbchain] = modelchain
			break
	    # Now re-write all of those ensemble files with the proper chains
	    # We also have to make sure that the chains appear in the same order in the ensemble as they
	    # do in the inputstruct (yes, lots of tedious constraints)
	    for pdbfile in pdblist:
		fout = open("temp", "w")
		for chain in inputchains:
		    fin = open(pdbfile, "r")
		    for aline in fin:
			if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and chainmap[aline[21]] == chain):
			    fout.write(aline[0:21] + chainmap[aline[21]] + aline[22:])
		    fin.close()
		fout.close()
		os.rename("temp", pdbfile)
    # First read the constraints file and see how tight the constraints are
    # Create a looser constraints file for round 1
    constraintsfile = "constraints.cst"
    try:
	fin = open(constraintsfile, "r")
	lowconstraintsfile = constraintsfile.split(".cst")[0] + "_low.cst"
	fout = open(lowconstraintsfile, "w")
	for aline in fin:
	    if (len(aline.strip()) == 0):
		continue
	    if ("BOUNDED" in aline):
		indx = aline.index("BOUNDED")
		lowbound = float(aline[indx:].split()[1])
		highbound = float(aline[indx:].split()[2])
		if (highbound < 10):
		    fout.write(aline[0:indx] + "BOUNDED " + str(lowbound) + " 12.0 " + aline[indx:].split()[3] + "\n")
		else:
		    fout.write(aline)
	    else:
		fout.write(aline)
	fin.close()
	fout.close()
    except:
	constraintsfile = None
    # Generate the flags file for round 1
    fout = open("flags", "w")
    fout.write("-in:path:database " + rosetta_db + "\n")
    fout.write("-in:file:s " + inputstruct + "\n")
    fout.write("-out:nstruct " + str(decoys) + "\n")
    if (randomize1):
	fout.write("-randomize1\n")
    if (randomize2):
	fout.write("-randomize2\n")
    if (ensemble1 or ensemble2):
	fout.write("-ensemble1 ensemble1\n")
	fout.write("-ensemble2 ensemble2\n")
    fout.write("-no_filters\n")
    fout.write("-score:weights talaris2013_cst\n")
    fout.write("-ignore_unrecognized_res\n")
    fout.write("-partners " + receptorchains + "_" + ligandchains + "\n")
    fout.write("-dock_pert 3 8 \n")
    if (constraintsfile):
	fout.write("-constraints:cst_file " + lowconstraintsfile + "\n")
    fout.write("-out:file:fullatom\n")
    fout.write("-overwrite\n")
    fout.write("-mute all\n")
    fout.close()
    # Now if there's an ensemble, we have to prepack the ensemble structures
    # We can use the same flags file, except change nstruct to 1
    if (ensemble1 or ensemble2):
	fin = open("flags", "r")
	fout = open("flags_prepack", "w")
	for aline in fin:
	    if (aline.startswith("-out:nstruct")):
		fout.write("-out:nstruct 1\n")
	    else:
		fout.write(aline)
	fin.close()
	fout.close()
	if (separateMPI_master):
	    commandline = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_prepack_protocol.mpi.linuxgccrelease @flags_prepack > dock.out) >& dock.err\""
	else:
	    commandline = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_prepack_protocol.mpi.linuxgccrelease @flags_prepack > dock.out) >& dock.err"
	print commandline
	res, output = commands.getstatusoutput(commandline)
	if (res):
	    print "ERROR: The Rosetta MPI prepack docker crashed!"
	    print "Output: " + output
	    exit()
	# Now we have to score those .ppk files and update the ensemble lists with the centroid and fa scores
	# Some versions of the prepacker do this automatically, so let's make sure Rosetta did what
	# it was supposed to do
	fin = open("ensemble1", "r")
	badRosetta = True
	for aline in fin:
	    try:
		float(aline)
		badRosetta = False
		break
	    except:
		pass
	fin.close()
	if (badRosetta):
	    commandline = "python ../../score_ensembles.py ensemble1 ensemble2"
	    print commandline
	    res, output = commands.getstatusoutput(commandline)
	    if (res):
		print "ERROR: The ensemble scorer crashed!"
		print "Output: " + output
		exit()
    # Do round 1 in MPI mode
    if (separateMPI_master):
	commandline = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_protocol.mpi.linuxgccrelease @flags > dock.out) >& dock.err\""
    else:
	commandline = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_protocol.mpi.linuxgccrelease @flags > dock.out) >& dock.err"
    print commandline
    res, output = commands.getstatusoutput(commandline)
    if (res):
	print "ERROR: The Rosetta MPI docker, round 1, crashed!"
	print "Output: " + output
	exit()
    # Move the files to the appropriate output folder
    pdbprefix = inputstruct.split(".pdb")[0]
    outputfiles = glob.glob(pdbprefix + "_*.pdb")
    if (outputdir != "."):
	try:
	    os.mkdir(outputdir)
	except:
	    pass
	for output in outputfiles:
	    os.rename(output, outputdir + "/" + output)
	outputfiles = glob.glob(outputdir + "/" + pdbprefix + "_*.pdb")
    # If there were no constraints, we're done
    if (constraintsfile):
	# Now let's find the top 100 PDBs by constraint-only score, using total_score to break ties (if no constraints are
	# given, then the ranking is basically just the total_score)
	totalS = []
	constraintS = []
	for outputfile in outputfiles:
	    fin = open(outputfile, "r")
	    cscore = 0.0
	    score = 9999999999.0
	    for aline in fin:
		if (aline.startswith("pose")):
		    score = float(aline.split()[len(aline.split())-1])
		elif (aline.startswith("atom_pair_constraint")):
		    cscore = float(aline.split()[len(aline.split())-1])
	    fin.close()
	    totalS.append(score)
	    constraintS.append(cscore)
	# Now sort
	for i in range(0, len(totalS)-1):
	    bestindx = i
	    for j in range(i, len(totalS)):
		if (constraintS[j] < constraintS[bestindx]):
		    bestindx = j
		elif (constraintS[j] == constraintS[bestindx] and totalS[j] < totalS[bestindx]):
		    bestindx = j
	    temp = outputfiles[i]
	    outputfiles[i] = outputfiles[bestindx]
	    outputfiles[bestindx] = temp
	    temp = totalS[i]
	    totalS[i] = totalS[bestindx]
	    totalS[bestindx] = temp
	    temp = constraintS[i]
	    constraintS[i] = constraintS[bestindx]
	    constraintS[bestindx] = temp
	winners = outputfiles[0:int(len(outputfiles) / 10)]
	# Generate a list of commands for doing all the round 2 docking simulations
	# Try to find the default, single processor Rosetta binary
	platform = "default.linuxgccrelease"
	binaries = glob.glob(rosetta_bin + "/docking_protocol.default.*")
	if (len(binaries) > 0):
	    platform = binaries[0].split("docking_protocol")[1]
	    platform = platform[1:] # Get rid of the leading .
	else:
	    binaries = glob.glob(rosetta_bin + "/docking_protocol.static.*")
	    if (len(binaries) > 0):
		platform = binaries[0].split("docking_protocol")[1]
		platform = platform[1:] # Get rid of the leading
	commandlines = []
	fout = open("commandslist", "w")
	for i in range(0, len(winners)):
	    # Generate the flags file
	    commandline = rosetta_bin + "/docking_protocol." + platform + " "
	    commandline = commandline + "-in:path:database " + rosetta_db + " "
	    commandline = commandline + "-in:file:s " + winners[i].strip() + " "
	    commandline = commandline + "-out:nstruct 10 "
	    commandline = commandline + "-no_filters "
	    if (ensemble1 or ensemble2):
		commandline = commandline + "-ensemble1 ensemble1 -ensemble2 ensemble2 "
	    commandline = commandline + "-score:weights talaris2013_cst "
	    commandline = commandline + "-ignore_unrecognized_res "
	    commandline = commandline + "-partners " + receptorchains + "_" + ligandchains + " "
	    commandline = commandline + "-dock_pert 3 8 "
	    commandline = commandline + "-constraints:cst_file " + constraintsfile + " " # Use the tight constraints now
	    commandline = commandline + "-out:file:fullatom "
	    commandline = commandline + "-overwrite"
	    fout.write(commandline + "\n")
	fout.close()
	# Spawn the child processes
	# For some reason, you cannot have mpi4py active in this script to spawn Rosetta processes
	# I think it has something to do with the fact that Rosetta is MPI-aware, so having this script running with
	# mpi4py activates MPI in the environment Rosetta starts in, so MPI-things get activated in Rosetta when they should
	# not be activated (since you are intending to only submit a lot of individual single-CPU Rosetta sessions) and
	# bad things happen
	# The method below seems to work
	#
	# Calculate what machines will be given to each machine, duplicating if a machine has multiple CPUs
	thishost = socket.gethostname()
	hosts = []
	fin = open(hostfile, "r")
	for aline in fin:
	    if (len(aline.strip()) == 0):
		continue
	    hostname = aline.split()[0]
	    ncpus = int(aline.split()[1].split("=")[1]) # hostname slots=ncpus max-slots=ncpus
	    for i in range(0, ncpus):
		hosts.append(hostname)
	fin.close()
	processes = []
	if (separateMPI_master):
	    commandline = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " python ../../dockR2.py > dock.out) >& dock.err\""
	else:
	    commandline = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " python ../../dockR2.py > dock.out) >& dock.err"
	res, output = commands.getstatusoutput(commandline)
	if (res):
	    print "ERROR: The Rosetta MPI docker, round 2, crashed!"
	    print "Output: " + output
	    exit()
	# Move the outputs if desired
	outputfiles = glob.glob(pdbprefix + "_*_*.pdb")
	if (outputdir != "."):
	    try:
		os.mkdir(outputdir + "/round2")
	    except:
		pass
	    pdbprefix = inputstruct.split(".pdb")[0]
	    for output in outputfiles:
		os.rename(output, outputdir + "/round2/" + output)
	    outputfiles = glob.glob(outputdir + "/round2/" + pdbprefix + "_*_*.pdb")
	else:
	    try:
		os.mkdir("round2")
	    except:
		pass
	    pdbprefix = inputstruct.split(".pdb")[0]
	    for output in outputfiles:
		os.rename(output, "round2/" + output)
	    outputfiles = glob.glob("round2/" + pdbprefix + "_*_*.pdb")
	# Now again, rank by increasing constraint-only score with total_score to break ties
	# Return the top 10 as the right answers
	totalS = []
	constraintS = []
    for outputfile in outputfiles:
	fin = open(outputfile, "r")
	cscore = 0.0
	score = 9999999999.0
	for aline in fin:
	    if (aline.startswith("pose")):
		score = float(aline.split()[len(aline.split())-1])
	    elif (aline.startswith("atom_pair_constraint")):
		cscore = float(aline.split()[len(aline.split())-1])
	fin.close()
	totalS.append(score)
	constraintS.append(cscore)
    # Now sort
    for i in range(0, len(totalS)-1):
	bestindx = i
	for j in range(i, len(totalS)):
	    if (constraintS[j] < constraintS[bestindx]):
		bestindx = j
	    elif (constraintS[j] == constraintS[bestindx] and totalS[j] < totalS[bestindx]):
		bestindx = j
	temp = outputfiles[i]
	outputfiles[i] = outputfiles[bestindx]
	outputfiles[bestindx] = temp
	temp = totalS[i]
	totalS[i] = totalS[bestindx]
	totalS[bestindx] = temp
	temp = constraintS[i]
	constraintS[i] = constraintS[bestindx]
	constraintS[bestindx] = temp
    winners = outputfiles[0:min(len(outputfiles), nfinal)]
    # Rename them so they are easy to find
    i = 1
    for winner in winners:
	if (outputdir != "."):
	    os.rename(winner, outputdir + "/final_%4.4i.pdb" % i)
	else:
	    os.rename(winner, "final_%4.4i.pdb" % i)
	i += 1

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
    if (len(aline.strip()) == 0 or aline.find("#") >= 0):
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
if (protocol == "msd"):
    addparams = True
    f = open("flags", "r")
    for aline in f:
	if ("-extra_res_fa" in aline):
	    addparams = False
	    break
    f.close()
else:
    addparams = False
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
try:
    os.remove(iRosetta_home + "/writing_to_queue")
except:
    pass
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
    # How many models do we need for certain protocols?
    if (protocol == "antibody"):
	fin = open("antibodyinput", "r")
	for aline in fin:
	    if (aline.startswith("NMODELS")):
		nmodels = int(aline.split("\t")[1])
		break
	fin.close()
    elif (protocol == "backrub"):
	fin = open("backrubinput", "r")
	for aline in fin:
	    if (aline.startswith("NMODELS")):
		nmodels = int(aline.split("\t")[1])
		break
	fin.close()
    elif (protocol == "dock"):
	fin = open("coarsedockinput", "r")
	for aline in fin:
	    if (aline.startswith("COARSEDECOYS")):
		nmodels = int(aline.split("\t")[1])
		break
	fin.close()
    elif (protocol == "pmutscan"):
	# Take as many as we can get
	nmodels = 99999
    # Query the backend nodes to find which nodes are free and submit to those nodes
    for host in hosts:
        hostname = host[0]
        slots = host[1]
        if (hostname == this_host):
	    process = Popen(args="uptime", stdout=PIPE, shell=True)
        elif (separateMPI_master):
            process = Popen(args="ssh " + MPI_master + " \'ssh " + hostname + " \"uptime\"\'", stdout=PIPE, shell=True)
        else:
            process = Popen(args="ssh " + hostname + " \"uptime\"", stdout=PIPE, shell=True)
        uptimeoutput = process.communicate()[0]
        if (hostname == this_host):
	    process = Popen(args="nproc", stdout=PIPE, shell=True)
        elif (separateMPI_master):
            process = Popen(args="ssh " + MPI_master + " \'ssh " + hostname + " \"nproc\"\'", stdout=PIPE, shell=True)
        else:
            process = Popen(args="ssh " + hostname + " \"nproc\"", stdout=PIPE, shell=True)
        nprocoutput = process.communicate()[0]
        try:
            load_1min = float(uptimeoutput.split()[len(uptimeoutput.split())-3].split(",")[0])
        except:
            load_1min = 0
	try:
            nproc = int(nprocoutput)
        except:
            nproc = 1
        if (float(load_1min) < 0.5 * float(nproc)):
            cpus_allocated = cpus_allocated + slots
            selected_hosts.append([hostname, slots])
	if (protocol in ["antibody", "dock", "pmutscan", "backrub"]):
	    if (cpus_allocated >= nmodels):
		cpus = nmodels
		runnable = True
		break
	    else:
		runnable = False
	else:
	    if (cpus_allocated >= cpus):
		runnable = True
		break
	    else:
		runnable = False
    os.remove(iRosetta_home + "/backend_query")
    if (protocol in ["antibody", "dock", "pmutscan", "backrub"]):
	# Ideally we want one processor per model, but if we have at least mincpus then let's go for it
	if (not(runnable) and cpus_allocated >= mincpus):
	    cpus = cpus_allocated
	    runnable = True
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
    if (protocol == "antibody"):
	# Generate the FASTA files
	fin = open("antibodyinput", "r")
	for aline in fin:
	    if (aline.startswith("LCSEQ")):
		fout = open("light.fasta", "w")
		fout.write("> Light Chain\n")
		fout.write(aline.split("\t")[1].strip() + "\n")
		fout.close()
	    elif (aline.startswith("HCSEQ")):
		fout = open("heavy.fasta", "w")
		fout.write("> Heavy Chain\n")
		fout.write(aline.split("\t")[1].strip() + "\n")
		fout.close()
	fin.close()
        if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../antibody_mpi.py MYAB " + str(nmodels) + " light.fasta heavy.fasta " + antibody_home + " " + rosetta_bin + " " + rosetta_db + " > antibody.out) >& antibody.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../antibody_mpi.py MYAB " + str(nmodels) + " light.fasta heavy.fasta " + antibody_home + " " + rosetta_bin + " " + rosetta_db + " > antibody.out) >& antibody.err"
    elif (protocol == "backrub"):
	# Generate the input PDB file and the weights file
	fin = open("backrubinput", "r")
	readingData = False
	for aline in fin:
	    if (aline[0:7] == "PDBFILE"):
		pdbfile = aline.split("\t")[1].strip()
		f2 = open(pdbfile, "w")
	    elif (aline[0:8] == "SCOREFXN"):
		weightsfile = aline.split("\t")[1].strip()
	    elif (aline[0:6] == "PARAMS"):
		paramsfile = aline.split("\t")[1].strip()
		f2 = open(paramsfile, "w")
	    elif (aline.startswith("BEGIN PDB")):
		readingData = True
	    elif (aline.startswith("END PDB")):
		f2.close()
		readingData = False
	    elif (aline.startswith("BEGIN PARAMS")):
		readingData = True
	    elif (aline.startswith("END PARAMS")):
		f2.close()
		readingData = False
	    elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
		weightsfile = "weights"
		f2 = open(weightsfile, "w")
		readingData = True
	    elif (aline[0:17] == "END SCOREFXN DATA"):
		f2.close()
		readingData = False
	    elif (readingData):
		f2.write(aline.strip() + "\n")
	    elif (aline[0:7] == "NMODELS"):
		ntemplates = int(aline.split("\t")[1].strip())
	    elif (aline[0:7] == "MAXRMSD"):
		maxRMSD = float(aline.split("\t")[1].strip())
	fin.close()
        if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../backrub.py " + pdbfile + " " + str(nmodels) + " " + str(maxRMSD) + " > backrub.out) >& backrub.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../backrub.py " + pdbfile + " " + str(nmodels) + " " + str(maxRMSD) + " > backrub.out) >& backrub.err"
    elif (protocol == "dock"):
	# Generate the input files
	fin = open("coarsedockinput", "r")
	readingData = False
	readingEnsemble = False
	for aline in fin:
	    if (aline[0:7] == "PDBFILE"):
		pdbfile = aline.split("\t")[1].strip()
		f2 = open(pdbfile, "w")
		if (readingEnsemble):
		    fensb.write(pdbfile.strip() + "\n")
	    elif (aline[0:8] == "SCOREFXN"):
		weightsfile = aline.split("\t")[1].strip()
	    elif (aline.startswith("BEGIN PDB")):
		readingData = True
		if (readingEnsemble):
		    f2 = open(aline.split()[2], "w")
		    fensb.write(aline.split()[2] + "\n")
	    elif (aline.startswith("END PDB")):
		f2.close()
		readingData = False
	    elif (aline.startswith("BEGIN ENSEMBLE")):
		readingEnsemble = True
		fensb = open("ensemble" + aline[14], "w")
	    elif (aline.startswith("END ENSEMBLE")):
		fensb.close()
		readingEnsemble = False
	    elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
		weightsfile = "weights"
		f2 = open(weightsfile, "w")
		readingData = True
	    elif (aline[0:17] == "END SCOREFXN DATA"):
		f2.close()
		readingData = False
	    elif (aline[0:14] == "BEGIN CST DATA"):
		f2 = open("constraints.cst", "w")
		readingData = True
	    elif (aline[0:12] == "END CST DATA"):
		f2.close()
		readingData = False
	    elif (readingData):
		f2.write(aline.strip() + "\n")
	    elif (aline[0:10] == "JUMPCONFIG"):
		jumpconfig = aline.split("\t")[1].strip()
	    elif (aline[0:6] == "ORIENT"):
		orient = aline.split("\t")[1].strip()
	    elif (aline[0:12] == "COARSEDECOYS"):
		ncoarse = int(aline.split("\t")[1].strip())
	    elif (aline[0:13] == "REFINEDDECOYS"):
		nrefined = int(aline.split("\t")[1].strip())
	fin.close()
	if (orient == "Global" or orient == "Fix Mov"):
	    randomize1 = True
	else:
	    randomize1 = False
	if (orient == "Global" or orient == "Fix Stat"):
	    randomize2 = True
	else:
	    randomize2 = False
	if (os.path.isfile("ensemble1")):
	    ensemble1 = True
	else:
	    ensemble1 = False
	if (os.path.isfile("ensemble2")):
	    ensemble2 = True
	else:
	    ensemble2 = False
	twoStageDock(cpus, pdbfile, "hostfile_" + jobID.strip(), jumpconfig.split("_")[0], jumpconfig.split("_")[1], ncoarse, nrefined, randomize1, randomize2, ensemble1, ensemble2)
    elif (protocol == "pmutscan"):
	# Get the pdbfile, resfile, and scorefxn from the input file
	fin = open("scaninput", "r")
	readingData = False
	for aline in fin:
	    if (aline[0:7] == "PDBFILE"):
		pdbfile = aline.split("\t")[1].strip()
		f2 = open(pdbfile, "w")
	    elif (aline[0:4] == "LIST"):
		listfile = aline.split("\t")[1].strip()
		f2 = open(listfile, "w")
	    elif (aline[0:8] == "SCOREFXN"):
		weightsfile = aline.split("\t")[1].strip()
	    elif (aline.startswith("MUTANTTYPE")):
		mutanttype = aline.split("\t")[1].strip()
	    elif (aline.startswith("BEGIN PDB")):
		readingData = True
	    elif (aline.startswith("END PDB")):
		f2.close()
		readingData = False
	    elif (aline.startswith("BEGIN LIST")):
		readingData = True
	    elif (aline.startswith("END LIST")):
		f2.close()
		readingData = False
	    elif (aline[0:19] == "BEGIN SCOREFXN DATA"):
		weightsfile = "weights"
		f2 = open(weightsfile, "w")
		readingData = True
	    elif (aline[0:17] == "END SCOREFXN DATA"):
		f2.close()
		readingData = False
	    elif (readingData):
		f2.write(aline.strip() + "\n")
	fin.close()
	# Write the flags file
	fout = open("flags", "w")
	fout.write("-in:file:s " + pdbfile + "\n")
	fout.write("-in:path:database " + rosetta_db + "\n")
	if (mutanttype == "DOUBLE"):
	    fout.write("-double_mutant_scan\n")
	fout.write("-mutants_list " + listfile + "\n")
	fout.write("-mute basic core\n")
	fout.close()
	if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/pmut_scan_parallel.mpi.linuxgccrelease @flags > pmut.out) >& pmut.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/pmut_scan_parallel.mpi.linuxgccrelease @flags > pmut.out) >& pmut.err"
    elif (protocol == "msd"):
	if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/mpi_msd.mpi.linuxgccrelease @flags > msd.out) >& msd.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/mpi_msd.mpi.linuxgccrelease @flags > msd.out) >& msd.err"
    if (protocol != "dock"):
	print command
	p = Popen(args=command, stdout=PIPE, shell=True)
	(out, err) = p.communicate()
# Now package the results up
if (protocol == "msd"):
    # Did we crash?  The answer is yes if no output files are present and msd.err has content in it
    crashed = False
    errmsg = ""
    outputfiles = glob.glob("msd_output*.pdb")
    if (len(outputfiles) == 0):
	crashed = True
    f = open("msd.err", "r")
    for aline in f:
        errmsg = errmsg + aline.strip() + "\n"
    f.close()
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
        exit()
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
elif (protocol == "antibody"):
    # Did we crash?  The answer is yes then we did not get the required number of output models
    crashed = False
    errmsg = ""
    f = open("antibody.err", "r")
    for aline in f:
        errmsg = errmsg + aline.strip() + "\n"
    f.close()
    # Sometimes warnings get outputted to the error file, so we can't know if an error really happening
    # simply by checking to see if antibody.err has content in it
    outputfiles = glob.glob("*.pdb")
    if (len(outputfiles) < nmodels):
	crashed = True
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
        exit()
elif (protocol == "dock"):
    # Did we crash?  The answer is yes then we did not get the required number of output models
    crashed = False
    errmsg = ""
    f = open("dock.err", "r")
    for aline in f:
        errmsg = errmsg + aline.strip() + "\n"
    f.close()
    # Sometimes warnings get outputted to the error file, so we can't know if an error really happening
    # simply by checking to see if dock.err has content in it
    outputfiles = glob.glob("final_????.pdb")
    if (len(outputfiles) == 0):
	crashed = True
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
        exit()
elif (protocol == "backrub"):
    # Did we crash?  The answer is yes then we did not get the required number of output models
    crashed = False
    errmsg = ""
    f = open("backrub.err", "r")
    for aline in f:
        errmsg = errmsg + aline.strip() + "\n"
    f.close()
    # Sometimes warnings get outputted to the error file, so we can't know if an error really happening
    # simply by checking to see if dock.err has content in it
    outputfiles = glob.glob(pdbfile.split(".pdb")[0] + "_*.pdb")
    if (len(outputfiles) == 0):
	crashed = True
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
        exit()
elif (protocol == "pmutscan"):
    # Did we crash?  The answer is yes if we have an error message in pmut.err
    crashed = False
    errmsg = ""
    f = open("pmut.err", "r")
    for aline in f:
        errmsg = errmsg + aline.strip() + "\n"
        if ("error" in aline.lower()):
	    crashed = True
    f.close()
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
    else:
	# Easy, just rename pmut.out so the client can find it
	fout = open("pmut.out", "a")
	fout.write("#MODEL\t" + pdbfile.strip() + "\n")
	fout.close()
	os.rename("pmut.out", "results.scan")
    exit()
archive = gzip.open("resultstemp.gz", "wb")
for pdbfile in outputfiles:
    archive.write("BEGIN PDB " + pdbfile.strip() + "\n")
    f = open(pdbfile.strip(), "r")
    for aline in f:
        archive.write(aline.strip() + "\n")
    f.close()
    archive.write("END PDB " + pdbfile.strip() + "\n")
archive.close()
if (protocol == "msd"):
    # MSD does not produce true ensembles, since it's going to have multiple states of the protein
    os.rename("resultstemp.gz", "results.msdar")
else:
    os.rename("resultstemp.gz", "results.ensb")
# Remove PDB files to keep the server tidy
pdbfiles = glob.glob("*.pdb")
for pdbfile in pdbfiles:
    os.remove(pdbfile)
try:
    shutil.rmtree("round2/" + ID, ignore_errors=True)
except:
    pass