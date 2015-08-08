### IMPORT WHATEVER YOU NEED HERE

### DO NOT CHANGE THE NAME OR LIST OF ARGUMENTS
### ANY DATA YOU NEED SHOULD BE IN THE INPUTFILE
def runJob(inputfile, 
	   protocolname,
	   jobID,
	   rosetta_bin,
	   rosetta_db,
	   iRosetta_home,
	   separateMPI_master,
	   MPI_master,
	   mpiexec,
	   hostfile_arg,
	   hostfile,
	   numproc_arg,
	   cpus):
    ### inputfile: The data uploaded in the GUI
    ### protocolname: The keyword for this protocol
    ### jobID: The ID given to this job when the user submitted it to this server
    ### rosetta_bin: The path to the Rosetta executables on this system
    ### rosetta_db: The path to the Rosetta database on this system
    ### iRosetta_home: The path to where this server runs the InteractiveROSETTA server (where
    ###                daemon_server lives)
    ### separateMPI_master: True if this computer is not the computer that runs parallel jobs,
    ###                     in this case the script sends the submission to the master via ssh
    ### MPI_master: The name of the computer that can run MPI if not the current one
    ### mpiexec: The command needed to run MPI on this system
    ### hostfile_arg: The hostfile argument keyword (e.g. -hostfile)
    ### hostfile: The file containing the list of hosts
    ### numproc_arg: The number of processors argument (e.g. -np)
    ### cpus: The number of nodes allocated for this job
    ### This script is running from <iRosetta_home>/results/<jobID>
    datadir = "../../modules/" + protocolname + "/data"
    # UNPACK THE INPUTFILE, ALL THE DATA NEEDED TO REGENERATE THE FILES LOCALLY SHOULD BE IN HERE
    fin = open(inputfile, "r")
    ### READ THE DATA, GENERATE THE FILES...
    ### ...
    ### ...
    fin.close()
    ### PREPARE THE INPUTS
    ### ...
    ### ...
    ### GENERATE THE FLAGS FILE MAYBE?
    ### ...
    ### ...
    ### EXAMPLE PROGRAM SUBMISSION:
    if (separateMPI_master):
	command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " " + hostfile + " " + numproc_arg + " 1 " + rosetta_bin + "/<prog> <args> > " + protocolname + ".out) >& " + protocolname + ".err\""
    else:
	command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " " + hostfile + " " + numproc_arg + " 1 " + rosetta_bin + "/<prog> <args> > " + protocolname + ".out) >& " + protocolname + ".err"
    print command
    p = Popen(args=command, stdout=PIPE, shell=True)
    (out, err) = p.communicate()
    ### YOU MAY NEED TO CALL MULTIPLE ROSETTA PROGRAMS AND THAT'S OKAY
    ### ...
    ### ...
    ### RETURN A LIST OF THE EXPECTED OUTPUTS SO THE CALLING SCRIPT CAN PACKAGE THEM
    ### THE CALLING SCRIPTS ASSUMES THAT THE FILES ARE IN <iRosetta_home>/results/<jobID>
    ### SO IF expected_outputs = ["myfile.pdb"], THE CALLING SCRIPT EXPECTS IT TO BE AT:
    ### <iRosetta_home>/results/<jobID>/myfile.pdb
    ### IF YOUR OUTPUT IS A LIST OF PDBS IN AN ENSEMBLE, USE "ensb" AS packagetype
    ### OTHERWISE, A GZ FILE WILL BE GENERATED WITH ALL YOUR FILES AND THIS WILL BE UNPACKED
    ### BY THE GUI
    expected_outputs = []
    ### ADD OUTPUTFILES
    ### ...
    ### ...
    packagetype = "ensb"
    return expected_outputs, packagetype