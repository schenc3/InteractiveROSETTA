import glob
import os
import os.path
import sys
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
iRosetta_home = "/var/www/html/walcob/InteractiveROSETTA/server"
# This is not the same file as for the Windows backend
# It is the list of nodes that are MPI-capable for running MPI C++ Rosetta
# mpi_hostlist is simply the hostlist file used for mpirun
mpi_hostlist = iRosetta_home + "/hostfile"
# This value is true if the server running daemon_server.py is not able to submit MPI jobs to a cluster
# If True, you need to specify the name of the machine that can submit MPI jobs so this script can ssh into it
separateMPI_master = False
MPI_master = "iR-testing"
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
# Location of the fragment tools package
fragment_tools = "/usr/local/rosetta/tools/fragment_tools"
# Location of the setup_RosettaCM.py script
cm_scripts = "/usr/local/rosetta/tools/protein_tools/scripts"
# ====================================================================================================================

# Attempt to import any custom modules
olddir = os.getcwd()
os.chdir("../../modules")
sys.path.append(os.getcwd())
modules = {}
for moduledir in glob.glob("*"):
    # Ignore the template
    if (moduledir == "template"):
        continue
    try:
        module = __import__(moduledir)
    except:
        print "Failed to import " + moduledir + ", does " + moduledir + "/__init__.py exist?"
        continue
    try:
        # Does the runJob function exist?
        module.runJob
    except:
        print "Failed to import " + moduledir + ", could not find the function \"runJob\""
        continue
    modules[moduledir] = module
os.chdir(olddir)

def AA3to1(resn):
    indx3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR HOH ".find(resn)
    if (indx3 < 0):
        return "Z"
    else:
        indx = indx3 / 4
        return "ACDEFGHIKLMNPQRSTVWYO"[indx]

def doDock(nproc, inputstruct, hostfile, receptorchains, ligandchains, decoys, nfinal, randomize1, randomize2, ensemble1, ensemble2):
    outputdir = "."
    writeStatus("Running")
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
            commandline = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_prepack_protocol.mpi.linuxgccrelease @flags_prepack > dock.out) &> dock.err\""
        else:
            commandline = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_prepack_protocol.mpi.linuxgccrelease @flags_prepack > dock.out) &> dock.err"
        print commandline
        res, output = commands.getstatusoutput(commandline)
        if (res):
            print "ERROR: The Rosetta MPI prepack docker crashed!"
            print "Output: " + output
            writeStatus("Failed")
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
                writeStatus("Failed")
                exit()
    # Do docking
    if (separateMPI_master):
        commandline = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_protocol.mpi.linuxgccrelease @flags > dock.out) &> dock.err\""
    else:
        commandline = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " -hostfile " + hostfile + " -np " + str(nproc) + " " + rosetta_bin + "/docking_protocol.mpi.linuxgccrelease @flags > dock.out) &> dock.err"
    print commandline
    res, output = commands.getstatusoutput(commandline)
    if (res):
        print "ERROR: The Rosetta MPI docker, round 1, crashed!"
        print "Output: " + output
        writeStatus("Failed")
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
    writeStatus("Finished")

def setupKIC():
    # Unpack everything we need from that inputfile and generate the flags file for KIC submission
    writeStatus("Running")
    fin = open("coarsekicinput", "r")
    nstruct = 1
    loopdata = []
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
        elif (aline[0:4] == "LOOP"):
            # Save the information in a list that we will iterate through later
            loopdata.append(aline.strip().split("\t")[1:])
        elif (aline[0:7] == "NSTRUCT"):
            nstruct = int(aline.split("\t")[1])
        elif (aline[0:7] == "PERTURB"):
            perturbType = aline.split("\t")[1].strip()
    fin.close()
    # Get the PDB data because we're probably going to have to modify it
    pdbdata = []
    fin = open(pdbfile, "r")
    last_res = "00000"
    for aline in fin:
        if (aline.startswith("ATOM") or aline.startswith("HETATM")):
            if (aline[22:27] != last_res):
                last_res = aline[22:27]
                pdbdata.append("")
            pdbdata[len(pdbdata)-1] += aline
    fin.close()
    # Read the dummy residues PDB data
    rsd_factory = []
    fin = open("../../data/residues.pdb", "r")
    for aline in fin:
        if (aline.startswith("ATOM") or aline.startswith("HETATM")):
            if (aline[22:27] != last_res):
                last_res = aline[22:27]
                rsd_factory.append("")
            rsd_factory[len(rsd_factory)-1] += aline
    fin.close()
    # Now construct dummy sequences for de novo loops
    for i in range(0, len(loopdata)):
        [loopType, sequence, begin, pivot, end] = loopdata[i]
        loopBegin = int(begin)
        pivot = int(pivot)
        loopEnd = int(end)
        if (loopType == "DE NOVO"):
            # Since this is a new sequence being added, we first have to delete all the residues
            # between the beginning and ending points
            oldlen = 0
            for ires in range(loopEnd-1, loopBegin, -1):
                pdbdata.pop(ires-1) # Rosetta indices are from 1, not 0
                #pose.delete_polymer_residue(ires)
                oldlen += 1
            # Now we have to add the sequence using our nifty little "rsd_factory" pose
            # The residues will have coordinates in weird places but it doesn't matter because
            # KIC fixes that and puts them in the right place; they don't need to start out anywhere
            # near being right
            offset = 0
            for AA in sequence.strip():
                indx = "ACDEFGHIKLMNPQRSTVWY".find(AA) + 1
                # Now we have to update the chain IDs for this new loop (currently blank)
                updatedstr = ""
                for entry in rsd_factory[indx].split("\n"):
                    updatedstr += entry[0:21] + pdbdata[loopBegin-1][21] + entry[22:] + "\n"
                pdbdata.insert(loopBegin+offset, updatedstr)
                #pose.append_polymer_residue_after_seqpos(Residue(rsd_factory.residue(indx)), loopBegin+offset, True)
                offset = offset + 1
            # Now we have to update the begin and end points of the other loops if necessary
            for j in range(0, len(loopdata)):
                if (loopBegin < int(loopdata[j][2])):
                    loopdata[j][2] = int(loopdata[j][2]) + len(sequence) - oldlen
                if (loopBegin < int(loopdata[j][3])):
                    loopdata[j][3] = int(loopdata[j][3]) + len(sequence) - oldlen
                if (loopBegin < int(loopdata[j][4])):
                    loopdata[j][4] = int(loopdata[j][4]) + len(sequence) - oldlen
            # Now maybe the sequence is longer than what was originally the length of the sequence
            # between start and end, so we need to recalculate the loop end
            loopEnd = loopBegin + len(sequence.strip()) + 1
        if (loopType == "DE NOVO"):
            # This has to be hard-coded, because the loop is not actually there until coarse modeling happens so there's no pivot point
            # other than the loop anchor residues
            loopdata[i][3] = loopEnd
    # Rewrite the PDB data
    fout = open(pdbfile, "w")
    for i in range(0, len(pdbdata)):
        data = pdbdata[i]
        # Renumber the residues
        for aline in data.split("\n"):
            if (len(aline.strip()) > 0):
                fout.write(aline[0:22] + "%4i" % (i+1) + aline[26:] + "\n")
    fout.close()
    # Write the loops file
    fout = open("in.loop", "w")
    for [loopType, sequence, begin, pivot, end] in loopdata:
        if (loopType == "REFINE"):
            fout.write("LOOP " + str(begin) + " " + str(end) + " " + str(pivot) + " 0 0\n")
        else:
            fout.write("LOOP " + str(begin) + " " + str(end) + " " + str(pivot) + " 0 1\n")
    fout.close()
    # Write the flags file
    fout = open("flags", "w")
    fout.write("-in:file:s " + pdbfile + "\n")
    fout.write("-database " + rosetta_db + "\n")
    fout.write("-loops:loop_file in.loop\n")
    fout.write("-loops:remodel perturb_kic\n")
    fout.write("-loops:refine refine_kic\n")
    fout.write("-ignore_unrecognized_res\n")
    fout.write("-overwrite\n")
    fout.write("-score:weights weights\n")
    fout.write("-out:nstruct " + str(nstruct) + "\n")
    fout.close()
    # Return the output stem so the ensemble packer can find the PDB outputs
    return pdbfile.split(".pdb")[0] + "_"
    writeStatus("Finished")

def doFlexPep():
    writeStatus("Running")
    # Generate the input files
    fin = open("flexpepinput", "r")
    readingData = False
    have_cst = False
    for aline in fin:
        if (aline[0:7] == "PDBFILE"):
            pdbfile = aline.split("\t")[1].strip()
            f2 = open(pdbfile, "w")
        elif (aline[0:8] == "SCOREFXN"):
            weightsfile = aline.split("\t")[1].strip()
        elif (aline.startswith("BEGIN PDB")):
            readingData = True
        elif (aline.startswith("END PDB")):
            f2.close()
            readingData = False
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
            have_cst = True
        elif (aline[0:12] == "END CST DATA"):
            f2.close()
            readingData = False
        elif (readingData):
            f2.write(aline.strip() + "\n")
        elif (aline.startswith("RECEPTOR")):
            receptorChains = aline.split("\t")[1].strip()
        elif (aline.startswith("PEPTIDE")):
            peptideChain = aline.split("\t")[1].strip()
        elif (aline.startswith("DECOYS")):
            ndecoys = int(aline.split("\t")[1].strip())
        elif (aline.startswith("RETURNMODELS")):
            nreturn = int(aline.split("\t")[1].strip())
        elif (aline.startswith("FLEXMODE")):
            mode = aline.split("\t")[1].strip()
    fin.close()
    if (mode == "ABINITIO"):
        # Generate the FASTA file for the peptide
        # This is needed so we can generate fragments
        fin = open(pdbfile, "r")
        fout = open("peptide.fasta", "w")
        fout.write("> Peptide Sequence\n")
        lastres = "0000"
        for aline in fin:
            if (aline.startswith("ATOM") and aline[21] == peptideChain and lastres != aline[22:26]):
                lastres = aline[22:26]
                fout.write(AA3to1(aline[17:20]))
        fout.write("\n")
        fin.close()
        fout.close()
        # Now we have to run PSIBLAST to get a position specific matrix to generate the sequence
        # profile for fragment selections
        commandline = "blastpgp -i peptide.fasta -j 2 -d nr.15 -C checkpoint.pssm"
        print commandline
        res, output = commands.getstatusoutput(commandline)
        if (res):
            print "ERROR: blastpgp failed!"
            print output
            writeStatus("Failed")
            exit()
        # Now we have to use Rosetta's script to convert the binary checkpoint.pssm to a text
        # file that can be used for sequence profiles
        commandline = "perl ../../gen-sequence-profile.pl peptide.fasta checkpoint.pssm frags.chk"
        print commandline
        res, output = commands.getstatusoutput(commandline)
        if (res):
            print "ERROR: gen-sequence-profile.pl failed!"
            print output
            writeStatus("Failed")
            exit()
        # Now run psipred to get the secondary structure prediction
        commandline = "runpsipred_single peptide.fasta"
        print commandline
        res, output = commands.getstatusoutput(commandline)
        if (res):
            print "ERROR: PSIPRED failed!"
            print output
        # Now generate a simple weighting scheme for fragment selection
        fout = open("simple.wghts", "w")
        fout.write("# score name        priority  wght   max_allowed  extras\n")
        fout.write("SecondarySimilarity 350 1.0 - psipred\n")
        fout.write("RamaScore 150 2.0 - psipred\n")
        fout.write("ProfileScoreL1 200 2.0 -\n")
        fout.write("SequenceIdentity 100 1.0 -") # No newline here otherwise you get a duplicate!!!
        fout.close()
        # Setup and run fragment generation
        fout = open("flags", "w")
        fout.write("-database " + rosetta_db + "\n")
        fout.write("-in::file::vall " + fragment_tools + "/vall.apr24.2008.extended.gz\n")
        fout.write("-in::file::fasta peptide.fasta\n")
        fout.write("-in::file::checkpoint frags.chk\n")
        fout.write("-frags::ss_pred peptide.ss2 psipred\n")
        fout.write("-frags::scoring::config simple.wghts\n")
        fout.write("-frags::bounded_protocol\n")
        fout.write("-frags::frag_sizes 3 5 9\n")
        fout.write("-frags::n_candidates 200\n")
        fout.write("-frags::n_frags 200\n")
        fout.write("-out::file::frag_prefix frags\n")
        fout.write("-frags::describe_fragments frags.fsc\n")
        fout.close()
        if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " 1 " + rosetta_bin + "/fragment_picker.mpi.linuxgccrelease @flags > frag.out) &> frag.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " 1 " + rosetta_bin + "/fragment_picker.mpi.linuxgccrelease @flags > frag.out) &> frag.err"
        print command
        p = Popen(args=command, stdout=PIPE, stderr=PIPE, shell=True)
        (out, err) = p.communicate()
        # As per the flexpepdock instructions, we have to shift the fragment positions in the
        # frags files by the number of residues in receptor chains
        # Count the number of receptor residues
        fin = open(pdbfile, "r")
        nResReceptor = 0
        last_res = "0000"
        for aline in fin:
            if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] in receptorChains and aline[22:26] != last_res):
                last_res = aline[22:26]
                nResReceptor += 1
        fin.close()
        fout = open("shift.csh", "w")
        fout.write("set ifragfile=\"frags.200.3mers\"\n")
        fout.write("set ofragfile=\"frags.200.shifted.3mers\"\n")
        fout.write("set nResReceptor=" + str(nResReceptor) + "\n")
        fout.write("awk '{if ( substr ( $0,1,3 ) == \"pos\" ) {print substr ( $0,0,18 ) sprintf (\"%4d\",substr ( $0,19,4 ) + '\"$nResReceptor\"' ) substr ( $0,23,1000 ) ; } else {print ; }}' $ifragfile > $ofragfile\n")
        fout.write("set ifragfile=\"frags.200.5mers\"\n")
        fout.write("set ofragfile=\"frags.200.shifted.5mers\"\n")
        fout.write("awk '{if ( substr ( $0,1,3 ) == \"pos\" ) {print substr ( $0,0,18 ) sprintf (\"%4d\",substr ( $0,19,4 ) + '\"$nResReceptor\"' ) substr ( $0,23,1000 ) ; } else {print ; }}' $ifragfile > $ofragfile\n")
        fout.write("set ifragfile=\"frags.200.9mers\"\n")
        fout.write("set ofragfile=\"frags.200.shifted.9mers\"\n")
        fout.write("awk '{if ( substr ( $0,1,3 ) == \"pos\" ) {print substr ( $0,0,18 ) sprintf (\"%4d\",substr ( $0,19,4 ) + '\"$nResReceptor\"' ) substr ( $0,23,1000 ) ; } else {print ; }}' $ifragfile > $ofragfile\n")
        fout.close()
        commandline = "tcsh shift.csh"
        print commandline
        res, output = commands.getstatusoutput(commandline)
        if (res):
            print "ERROR: Could not shift the frag files!"
            print output
            writeStatus("Failed")
            exit()
    # Now that we have the fragments, we can run flexpepdock
    # Generate the new frags file
    fout = open("flags", "w")
    fout.write("-database " + rosetta_db + "\n")
    fout.write("-s " + pdbfile + "\n")
    fout.write("-ex1\n")
    fout.write("-ex2aro\n")
    fout.write("-use_input_sc\n")
    if (mode == "ABINITIO"):
        fout.write("-frag3 frags.200.shifted.3mers\n")
        fout.write("-frag5 frags.200.shifted.5mers\n")
        fout.write("-frag9 frags.200.shifted.9mers\n")
        fout.write("-lowres_abinitio\n")
    fout.write("-pep_refine\n")
    fout.write("-nstruct " + str(ndecoys) + "\n")
    if (have_cst):
        fout.write("-constraints:cst_file constraints.cst\n")
        if (mode == "ABINITIO"):
            fout.write("-place_peptide_on_binding_site\n")
    fout.write("-receptor_chain " + receptorChains + "\n")
    fout.write("-peptide_chain " + peptideChain + "\n")
    fout.write("-score:weights talaris2013_cst\n")
    fout.close()
    # Run it
    if (separateMPI_master):
        command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/FlexPepDocking.mpi.linuxgccrelease @flags > flexpep.out) &> flexpep.err\""
    else:
        command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/FlexPepDocking.mpi.linuxgccrelease @flags > flexpep.out) &> flexpep.err"
    print command
    p = Popen(args=command, stdout=PIPE, stderr=PIPE, shell=True)
    (out, err) = p.communicate()
    # Okay, now let's take the top N models as ranked by the constraints score first and
    # reweighted FlexPepDock score second
    # We can get this information out of the score file
    pdbstem = pdbfile.split(".pdb")[0]
    pdblist = []
    cst_scores = []
    rw_scores = []
    for i in range(0, ndecoys):
        pdblist.append(pdbstem + "_%4.4i.pdb" % (i+1))
        cst_scores.append(0.0)
        rw_scores.append(0.0)
    fin = open("score.sc", "r")
    for aline in fin:
        if ("total_score" in aline):
            fields = aline.split()
            try:
                rw_indx = fields.index("reweighted_sc")
            except:
                rw_indx = fields.index("total_score")
            try:
                cst_indx = fields.index("atom_pair_constraint")
            except:
                cst_indx = -1
            desc_indx = fields.index("description")
        elif (aline.startswith("SCORE:")):
            try:
                pdb = aline.split()[desc_indx] + ".pdb"
                score = float(aline.split()[rw_indx])
                rw_scores[pdblist.index(pdb)] = score
                if (cst_indx >= 0):
                    score = float(aline.split()[cst_indx])
                    cst_scores[pdblist.index(pdb)] = score
            except:
                pass
    fin.close()
    # Now do the select sorting
    for i in range(0, len(pdblist)-1):
        best = i
        for j in range(i+1, len(pdblist)):
            if (cst_scores[j] < cst_scores[best]):
                best = j
            elif (cst_scores[j] == cst_scores[best] and rw_scores[j] < rw_scores[best]):
                best = j
        temp = pdblist[best]
        pdblist[best] = pdblist[i]
        pdblist[i] = temp
        temp = rw_scores[best]
        rw_scores[best] = rw_scores[i]
        rw_scores[i] = temp
        temp = cst_scores[best]
        cst_scores[best] = cst_scores[i]
        cst_scores[i] = temp
    # Copy them to new filename
    for i in range(0, nreturn):
        outfile = "final_flexpep_%4.4i.pdb" % (i+1)
        os.system("cp " + pdblist[i] + " " + outfile)
    writeStatus("Finished")

def writeStatus(status):
    statusOut = open("status",'w+')
    statusOut.write(status)
    statusOut.close()

def doThread(hostfile):
    # Unpack the input files
    fin = open("threadinput", "r")
    writeStatus("Running")
    readingData = False
    templates = []
    targaligns = []
    tempaligns = []
    for aline in fin:
        if (aline.startswith("FASTA")):
            FASTAseq = aline.split("\t")[1].strip()
        elif (aline.startswith("PDBFILE")):
            pdbfile = aline.split("\t")[1].strip()
            templates.append(pdbfile)
            fout = open(pdbfile, "w")
        elif (aline.startswith("BEGIN PDB DATA")):
            readingData = True
        elif (aline.startswith("END PDB DATA")):
            fout.close()
            readingData = False
        elif (aline.startswith("TARGALIGN")):
            targaligns.append(aline.split("\t")[1].strip())
        elif (aline.startswith("TEMPALIGN")):
            tempaligns.append(aline.split("\t")[1].strip())
        elif (aline.startswith("NMODELS")):
            nstruct = int(aline.split("\t")[1])
        elif (readingData):
            fout.write(aline)
    fin.close()
    # Let's get the FASTA file
    fout = open("target.fasta", "w")
    fout.write("> target\n" + FASTAseq.strip() + "\n")
    fout.close()
    # First we have to generate fragments for the comparative modeler
    # Now we have to run PSIBLAST to get a position specific matrix to generate the sequence
    # profile for fragment selections
    commandline = "blastpgp -i target.fasta -j 2 -d nr.15 -C checkpoint.pssm"
    print commandline
    res, output = commands.getstatusoutput(commandline)
    if (res):
        print "ERROR: blastpgp failed!"
        print output
        writeStatus("Failed")
        exit()
    # Now we have to use Rosetta's script to convert the binary checkpoint.pssm to a text
    # file that can be used for sequence profiles
    commandline = "perl ../../gen-sequence-profile.pl target.fasta checkpoint.pssm frags.chk"
    print commandline
    res, output = commands.getstatusoutput(commandline)
    if (res):
        print "ERROR: gen-sequence-profile.pl failed!"
        print output
        writeStatus("Failed")
        exit()
    # Now run psipred to get the secondary structure prediction
    commandline = "runpsipred_single target.fasta"
    print commandline
    res, output = commands.getstatusoutput(commandline)
    if (res):
        print "ERROR: PSIPRED failed!"
        print output
        writeStatus("Failed")
        exit()
    # Now generate a simple weighting scheme for fragment selection
    fout = open("simple.wghts", "w")
    fout.write("# score name        priority  wght   max_allowed  extras\n")
    fout.write("SecondarySimilarity 350 1.0 - psipred\n")
    fout.write("RamaScore 150 2.0 - psipred\n")
    fout.write("ProfileScoreL1 200 2.0 -\n")
    fout.write("SequenceIdentity 100 1.0 -") # No newline here otherwise you get a duplicate!!!
    fout.close()
    # Setup and run fragment generation
    fout = open("flags", "w")
    fout.write("-database " + rosetta_db + "\n")
    fout.write("-in::file::vall " + fragment_tools + "/vall.apr24.2008.extended.gz\n")
    fout.write("-in::file::fasta target.fasta\n")
    fout.write("-in::file::checkpoint frags.chk\n")
    fout.write("-frags::ss_pred target.ss2 psipred\n")
    fout.write("-frags::scoring::config simple.wghts\n")
    fout.write("-frags::bounded_protocol\n")
    fout.write("-frags::frag_sizes 3 5 9\n")
    fout.write("-frags::n_candidates 200\n")
    fout.write("-frags::n_frags 200\n")
    fout.write("-out::file::frag_prefix frags\n")
    fout.write("-frags::describe_fragments frags.fsc\n")
    fout.close()
    if (separateMPI_master):
        command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " " + hostfile + " " + numproc_arg + " 1 " + rosetta_bin + "/fragment_picker.mpi.linuxgccrelease @flags > frag.out) &> frag.err\""
    else:
        command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " " + hostfile + " " + numproc_arg + " 1 " + rosetta_bin + "/fragment_picker.mpi.linuxgccrelease @flags > frag.out) &> frag.err"
    print command
    p = Popen(args=command, stdout=PIPE,stderr=PIPE, shell=True)
    (out, err) = p.communicate()
    # Okay, now we need to set up a multiple sequence alignment for the comparative modeler
    # The alignments as they are now are correct, but the target sequences may have gaps in
    # different places so all we have to do is line up all the target sequences and keep the
    # template sequences aligned to the same residues

    indx = 0
    while (True):
        getout = True
        for i in range(0, len(targaligns)):
            if (indx < len(targaligns[i])):
                getout = False
                break
        if (getout):
            break
        # Do we have a gap at this index in any of the target alignments?
        havegap = False
        for i in range(0, len(targaligns)):
            if (indx >= len(targaligns[i]) or targaligns[i][indx] == "-"):
                havegap = True
        if (havegap):
            # Yes, we do.  Add one into each target sequence that does not have a gap
            # and also add a gap to their partner alignments
            for i in range(0, len(targaligns)):
                if (indx >= len(targaligns[i])):
                    targaligns[i] += "-"
                    tempaligns[i] += "-"
                elif (targaligns[i][indx] != "-"):
                    targaligns[i] = targaligns[i][0:indx] + "-" + targaligns[i][indx:]
                    tempaligns[i] = tempaligns[i][0:indx] + "-" + tempaligns[i][indx:]
        indx += 1
    # Write the alignment out as a FASTA file
    fout = open("comparative.aln", "w")
    fout.write("> target\n" + targaligns[0].strip() + "\n\n")
    for i in range(0, len(tempaligns)):
        fout.write("> " + templates[i].split(".pdb")[0] + "\n" + tempaligns[i].strip() + "\n\n")
    fout.close()
    # Okay, now let's run the setup_RosettaCM.py script to generate some input files for CM
    commandline = "python " + cm_scripts + "/setup_RosettaCM.py --fasta target.fasta " + \
        "--alignment comparative.aln --alignment_format fasta --rosetta_bin " + rosetta_bin + " --templates "
    for pdbfile in templates:
        commandline += pdbfile.strip() + " "
    print commandline
    # For some reason the script below sometimes does not generate the appropriate pdbfiles
    # Let's help it along now
    try:
        os.mkdir("rosetta_cm")
    except:
        pass
    for template in templates:
        shutil.copy(template, "rosetta_cm/" + template.split(".pdb")[0] + "_201.pdb") # This seems to be what it looks for
    res, output = commands.getstatusoutput(commandline)
    if (res):
        print "%i ERROR: setup_RosettaCM.py failed!" %(res)
        print output
        print
        writeStatus("Failed")
    # This generates a folder called "rosetta_cm" that contains most of our input files
    os.chdir("rosetta_cm")
    # The flags file does not have the frags files in it, so let's add them
    fout = open("flags", "a")
    fout.write("-frag3 ../frags.200.3mers\n")
    fout.write("-frag5 ../frags.200.5mers\n")
    fout.write("-frag9 ../frags.200.9mers\n")
    fout.close()
    # There's also a score function issue, apparently cart_bonded and pro_close cannot be used
    # together, so turn off pro_close (https://www.rosettacommons.org/content/multiple-homology-modelling, pose #28)
    data = []
    fin = open("stage3_rlx.wts", "r")
    for aline in fin:
        if ("pro_close" not in aline):
            data.append(aline)
    fin.close()
    fout = open("stage3_rlx.wts", "w")
    for aline in data:
        fout.write(aline)
    fout.close()
    # Now we can run rosettaCM
    if (separateMPI_master):
        command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "/rosetta_cm; (" + mpiexec + " " + hostfile_arg + " ../hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/rosetta_scripts.mpi.linuxgccrelease @flags -database " + rosetta_db + " -nstruct " + str(nstruct) + " > cm.out) &> cm.err\""
    else:
        command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "/rosetta_cm; (" + mpiexec + " " + hostfile_arg + " ../hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/rosetta_scripts.mpi.linuxgccrelease @flags -database " + rosetta_db + " -nstruct " + str(nstruct) + " > cm.out) &> cm.err"
    print command
    p = Popen(args=command, stdout=PIPE,stderr=PIPE, shell=True)
    (out, err) = p.communicate()

    # Move the output files back out of this directory
    expected_outputs = []
    for i in range(0, nstruct):
        outputname = "final_cm_%4.4i.pdb" % (i+1)
        #os.rename("S_%4.4i.pdb" % (i+1), "../" + outputname)
        status,output = commands.getstatusoutput("cp -v S_%4.4i.pdb ../%s"%(i+1,outputname))
        print status
        print output
    os.chdir("..")
    print "doThread finished!"
    writeStatus("Finished")

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
    paramsfiles = glob.glob("*.fa.params")
    if (len(paramsfiles) > 0):
        f = open("flags", "a")
        f.write("-extra_res_fa")
        for params in paramsfiles:
            f.write(" " + params.strip())
        f.write("\n")
        f.close()
# Enter the queue
while (True):
    while (os.path.isfile(iRosetta_home + "/writing_to_queue")):
        time.sleep(10)
    f = open(iRosetta_home + "/writing_to_queue", "w")
    f.write(this_host)
    f.close()
    # Wait a while, in case other nodes got in as well, then check the hostname in the file
    # and wait some more if it is not this hostname
    time.sleep(3)
    fin = open(iRosetta_home + "/writing_to_queue", "r")
    filehostname = fin.readlines()[0].strip()
    fin.close()
    if (filehostname == this_host):
        break
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
    while (True):
        while (os.path.isfile(iRosetta_home + "/writing_to_queue")):
            time.sleep(10)
        f = open(iRosetta_home + "/writing_to_queue", "w")
        f.write(this_host)
        f.close()
        # Wait a while, in case other nodes got in as well, then check the hostname in the file
        # and wait some more if it is not this hostname
        time.sleep(3)
        fin = open(iRosetta_home + "/writing_to_queue", "r")
        filehostname = fin.readlines()[0].strip()
        fin.close()
        if (filehostname == this_host):
            break
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
    else:
        # Is the jobfile still there?  If not, get rid of this entry so we can proceed
        if (len(glob.glob(iRosetta_home + "/jobfiles/*" + data[0].strip() + "*")) == 0):
            f = open(iRosetta_home + "/rosetta_queue", "w")
            for aline in data[1:]:
                f.write(aline.strip() + "\n")
            f.close()
    os.remove(iRosetta_home + "/writing_to_queue")
    time.sleep(10)
# Loop for finding free nodes
while (True):
    while (True):
        while (os.path.isfile(iRosetta_home + "/backend_query")):
            time.sleep(10)
        f = open(iRosetta_home + "/backend_query", "w")
        f.write(this_host)
        f.close()
        # Wait a while, in case other nodes got in as well, then check the hostname in the file
        # and wait some more if it is not this hostname
        time.sleep(3)
        fin = open(iRosetta_home + "/backend_query", "r")
        filehostname = fin.readlines()[0].strip()
        fin.close()
        if (filehostname == this_host):
            break
    # How many models do we need for certain protocols?
    forced_minimum = False
    custom_module = False
    if (protocol == "msd"):
        nmodels = cpus
    elif (protocol == "antibody"):
        fin = open("antibodyinput", "r")
        for aline in fin:
            if (aline.startswith("NMODELS")):
                nmodels = int(aline.split("\t")[1])
                break
        fin.close()
    elif (protocol == "thread"):
        fin = open("threadinput", "r")
        for aline in fin:
            if (aline.startswith("NMODELS")):
                nmodels = int(aline.split("\t")[1])
                break
        fin.close()
    elif (protocol == "kic"):
        fin = open("coarsekicinput", "r")
        for aline in fin:
            if (aline.startswith("NSTRUCT")):
                nmodels = int(aline.split("\t")[1]) + 1
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
    elif (protocol == "flexpep"):
        fin = open("flexpepinput", "r")
        for aline in fin:
            if (aline.startswith("DECOYS")):
                nmodels = int(aline.split("\t")[1])
                break
        fin.close()
        # We might want to change this later...It's set up to get enough nodes that FlexPepDock
        # should not take more than a day to complete so it doesn't hold other people's jobs up
        # It shouldn't be an issue unless the user requested a bazillion models
        required_nodes = nmodels / 300
        availNodes = 0
        for host in hosts:
            availNodes += int(host[1])
        if (availNodes < required_nodes):
            fout = open("errreport", "w")
            fout.write("We don't have enough CPUs to complete your job, so we are aborting.\n")
            fout.close()
            os.remove(iRosetta_home + "/backend_query")
            exit()
    elif (protocol == "pmutscan"):
        # Take as many as we can get
        nmodels = 99999
    else:
        # Custom module
        custom_module = True
        installed = False
        # Is this protocol installed?
        for installedmodule in glob.glob("../../modules/*"):
            if (installedmodule.endswith(protocol)):
                installed = True
                break
        if (not(installed)):
            fout = open("errreport", "w")
            fout.write("The protocol " + protocol + " is not installed on this server, so we are aborting.\n")
            fout.close()
            os.remove(iRosetta_home + "/backend_query")
            exit()
        nmodels = -1
        fin = open(protocol + "input", "r")
        for aline in fin:
            if (aline.startswith("REQUESTED_NODES")):
                nmodels = int(aline.split("\t")[1])
            elif (aline.startswith("REQUIRED_NODES")):
                nmodels = int(aline.split("\t")[1])
                forced_minimum = True
                break
        fin.close()
        # Error if the GUI didn't tell us how many nodes to use
        if (nmodels <= 0):
            fout = open("errreport", "w")
            fout.write("Your protocol did not specify how many nodes to use, so we aborted the job.\n")
            fout.close()
            os.remove(iRosetta_home + "/backend_query")
            exit()
        # See how many nodes are available total and make sure we can get that number
        availNodes = 0
        for host in hosts:
            availNodes += int(host[1])
        if (forced_minimum and availNodes < nmodels):
            fout = open("errreport", "w")
            fout.write("We don't have enough CPUs to complete your job, so we are aborting.\n")
            fout.close()
            os.remove(iRosetta_home + "/backend_query")
            exit()
    # Query the backend nodes to find which nodes are free and submit to those nodes
    for host in hosts:
        # Wait a minute in case a new job is starting, we need to wait to let the CPUs
        # get busy enough that uptime reports accurate loads
        #time.sleep(60)
        hostname = host[0]
        slots = host[1]
        if (hostname == this_host):
            process = Popen(args="uptime", stdout=PIPE,stderr = PIPE, shell=True)
        elif (separateMPI_master):
            process = Popen(args="ssh " + MPI_master + " \'ssh " + hostname + " \"uptime\"\'", stdout=PIPE, stderr = PIPE, shell=True)
        else:
            process = Popen(args="ssh " + hostname + " \"uptime\"", stdout=PIPE, stderr = PIPE, shell=True)
        uptimeoutput = process.communicate()[0]
        if (hostname == this_host):
            process = Popen(args="nproc", stdout=PIPE, stderr = PIPE, shell=True)
        elif (separateMPI_master):
            process = Popen(args="ssh " + MPI_master + " \'ssh " + hostname + " \"nproc\"\'", stdout=PIPE, stderr = PIPE, shell=True)
        else:
            process = Popen(args="ssh " + hostname + " \"nproc\"", stdout=PIPE, stderr = PIPE, shell=True)
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
        if (cpus_allocated >= nmodels):
            cpus = nmodels
            runnable = True
            break
        else:
            runnable = False
    os.remove(iRosetta_home + "/backend_query")
    if (protocol in ["antibody", "dock", "pmutscan", "backrub", "kic", "thread"] or (custom_module and not(forced_minimum))):
        # Ideally we want one processor per model, but if we have at least mincpus then let's go for it
        if (not(runnable) and cpus_allocated >= mincpus):
            cpus = cpus_allocated
            runnable = True
    elif (protocol == "flexpep"):
        # This one really needs a lot to finish in time, this is calculated such that it should
        # get enough nodes to finish in one day at the most
        if (not(runnable) and cpus_allocated >= nmodels / 300):
            cpus = cpus_allocated
            runnable = True
    if (runnable):
        break
    selected_hosts = []
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
    packagetype = "default"
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
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../antibody_mpi.py MYAB " + str(nmodels) + " light.fasta heavy.fasta " + antibody_home + " " + rosetta_bin + " " + rosetta_db + " > antibody.out) &> antibody.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../antibody_mpi.py MYAB " + str(nmodels) + " light.fasta heavy.fasta " + antibody_home + " " + rosetta_bin + " " + rosetta_db + " > antibody.out) &> antibody.err"
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
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../backrub.py " + pdbfile + " " + str(nmodels) + " " + str(maxRMSD) + " > backrub.out) &> backrub.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " python ../../backrub.py " + pdbfile + " " + str(nmodels) + " " + str(maxRMSD) + " > backrub.out) &> backrub.err"
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
        doDock(cpus, pdbfile, "hostfile_" + jobID.strip(), jumpconfig.split("_")[0], jumpconfig.split("_")[1], ncoarse, nrefined, randomize1, randomize2, ensemble1, ensemble2)
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
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/pmut_scan_parallel.mpi.linuxgccrelease @flags > pmut.out) &> pmut.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/pmut_scan_parallel.mpi.linuxgccrelease @flags > pmut.out) &> pmut.err"
    elif (protocol == "kic"):
        outputstem = setupKIC()
        if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/loopmodel.mpi.linuxgccrelease @flags > kic.out) &> kic.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/loopmodel.mpi.linuxgccrelease @flags > kic.out) &> kic.err"
    elif (protocol == "msd"):
        if (separateMPI_master):
            command = "ssh " + MPI_master + " \"cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/mpi_msd.mpi.linuxgccrelease @flags > msd.out) &> msd.err\""
        else:
            command = "cd " + iRosetta_home + "/results/" + jobID.strip() + "; (" + mpiexec + " " + hostfile_arg + " hostfile_" + jobID.strip() + " " + numproc_arg + " " + str(cpus) + " " + rosetta_bin + "/mpi_msd.mpi.linuxgccrelease @flags > msd.out) &> msd.err"
    elif (protocol == "flexpep"):
        doFlexPep()
    elif (protocol == "thread"):
        doThread("hostfile_" + jobID.strip())
    else:
        try:
            (expected_outputs, packagetype) = modules[protocol].runJob(protocol + "input",
                                                 protocol,
                                                 jobID.strip(),
                                                 rosetta_bin,
                                                 rosetta_db,
                                                 iRosetta_home,
                                                 separateMPI_master,
                                                 MPI_master,
                                                 mpiexec,
                                                 hostfile_arg,
                                                 "hostfile_" + jobID.strip(),
                                                 numproc_arg,
                                                 cpus)
        except Exception as e:
            print "rosetta_submit.py - Cannot run protocol " + protocol
            fout = open("errreport", "w")
            fout.write("Could not run protocol " + protocol + "\n")
            fout.write(e.message)
            fout.close()
            exit()
    if (protocol not in ["dock", "flexpep", "thread"] and not(custom_module)):
        print command
        p = Popen(args=command, stdout=PIPE, stderr = PIPE, shell=True)
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
elif (protocol == "flexpep"):
    # Did we crash?  The answer is yes then we did not get the required number of output models
    crashed = False
    errmsg = ""
    try:
        f = open("flexpep.err", "r")
        for aline in f:
            errmsg = errmsg + aline.strip() + "\n"
        f.close()
    except:
        f = open("frag.err", "r")
        for aline in f:
            errmsg = errmsg + aline.strip() + "\n"
        f.close()
    # Sometimes warnings get outputted to the error file, so we can't know if an error really happening
    # simply by checking to see if dock.err has content in it
    outputfiles = glob.glob("final_flexpep_????.pdb")
    if (len(outputfiles) == 0):
        crashed = True
    if (crashed):
        f = open("errreport", "w")
        f.write(errmsg)
        f.close()
        exit()
elif (protocol == "thread"):
    # Did we crash?  The answer is yes then we did not get the required number of output models
    crashed = False
    errmsg = ""
    try:
        f = open("cm.err", "r")
        for aline in f:
            errmsg = errmsg + aline.strip() + "\n"
        f.close()
    except:
        f = open("frag.err", "r")
        for aline in f:
            errmsg = errmsg + aline.strip() + "\n"
        f.close()
    # Sometimes warnings get outputted to the error file, so we can't know if an error really happening
    # simply by checking to see if dock.err has content in it
    outputfiles = glob.glob("final_cm_????.pdb")
    print "outputfiles"
    print outputfiles
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
elif (protocol == "kic"):
    # Did we crash?  The answer is yes then we did not get the required number of output models
    crashed = False
    errmsg = ""
    f = open("kic.err", "r")
    for aline in f:
        errmsg = errmsg + aline.strip() + "\n"
    f.close()
    # Sometimes warnings get outputted to the error file, so we can't know if an error really happening
    # simply by checking to see if dock.err has content in it
    outputfiles = glob.glob(outputstem + "*.pdb")
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
else:
    # Custom module
    # The runJob argument should return a list of the expected output files
    # If it crashed, it should have already written the error message and exited by now
    outputfiles = expected_outputs
archive = gzip.open("resultstemp.gz", "wb")
for pdbfile in outputfiles:
    if (packagetype.lower() != "default" and packagetype.lower() != "ensb"):
        archive.write("BEGIN FILE " + pdbfile.strip() + "\n")
        f = open(pdbfile.strip(), "r")
        for aline in f:
            archive.write(aline.strip() + "\n")
        f.close()
        archive.write("END FILE " + pdbfile.strip() + "\n")
    else:
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
elif (packagetype.lower() == "default" or packagetype.lower() == "ensb"):
    os.rename("resultstemp.gz", "results.ensb")
else:
    os.rename("resultstemp.gz", "results.gz")
# Remove PDB files to keep the server tidy
pdbfiles = glob.glob("*.pdb")
for pdbfile in pdbfiles:
    os.remove(pdbfile)
try:
    shutil.rmtree("round2/" + ID, ignore_errors=True)
except:
    pass
