from mpi4py import MPI
import sys
import os
import glob
import subprocess
import shutil

def cleanPDB(pdbfile):
    # This function will look for and remove duplicate atoms
    # It will permanently modify the PDB that the user loaded
    # This shouldn't cause problems for other programs though
    data = []
    taken_nums = {}
    num_mapping = {}
    # Sometimes there are PDBs that have multiple residues with the same residue index
    # BioPython drops these, but it can lead to all kinds of problems later on
    # So I will keep a record of the backbone atoms of the current residue and if we encounter
    # a BB atom with the same residue index, we will assume it's a new residue and renumber the residue
    lastBBatoms = []
    altlocs_taken = ""
    f = open(pdbfile.strip(), "r")
    curr_res = "   0"
    for aline in f:
        if (aline.startswith("ATOM") or aline.startswith("HETATM")):
            res = aline[22:27] # Includes residue indx + the optional alternate letter
            if (res[0:4] != curr_res[0:4]): # New residue indx
                altlocs_taken = res[4] # Reset the taken altlocs
                curr_res = res
                atomtypes = []
                lastBBatoms = []
            # This is only done if this is a new residue, but not a new residue indx
            if (aline[22:27] != curr_res or aline[12:16] in lastBBatoms):
                curr_res = res
                atomtypes = []
                lastBBatoms = []
                # Assign the altloc to whatever the most recent altloc used was
                for char in " ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                    if (not(char in altlocs_taken)):
                        altlocs_taken = char + altlocs_taken
                        break
            res = res[0:4] + altlocs_taken[0]
            atomtype = aline[12:16]
            if (atomtype in [" C  ", " CA ", " O  ", " N  "]):
                lastBBatoms.append(atomtype)
            if (atomtype in atomtypes):
                # Find a new type for this atom
                stem = atomtype[0:2]
                for i in range(1, 100):
                    if (i < 10):
                        newtype = stem + str(i) + " "
                    else:
                        newtype = stem + str(i)
                    if (not(newtype in atomtypes)):
                        atomtypes.append(newtype)
                        break
                aline = aline[0:12] + newtype + aline[16:]
            else:
                atomtypes.append(atomtype)
            # Now check for alternate forms of residues (i.e. 30A, 30B, etc.)
            # Rename these so the each have a unique number, the user can delete extras later
            chain = aline[21]
            if (not(chain in taken_nums.keys())):
                taken_nums[chain] = []
            if (not(res[0:4] in taken_nums[chain])):
                taken_nums[chain].append(res[0:4])
                num_mapping[chain+res] = res[0:4]
            else:
                try:
                    aline = aline[0:22] + num_mapping[chain+res] + " " + aline[27:]
                except:
                    # Find a new ID
                    lastnum = int(taken_nums[chain][len(taken_nums[chain])-1]) + 1
                    num_mapping[chain+res] = "%4i" % lastnum
                    taken_nums[chain].append("%4i" % lastnum)
                    aline = aline[0:22] + num_mapping[chain+res] + " " + aline[27:]
            #data.append(aline.strip())
        #else:
        data.append(aline.strip())
    f.close()
    # Now write the updated data out
    f = open(pdbfile.strip(), "w")
    for aline in data:
        f.write(aline + "\n")
    f.close()

if (len(sys.argv) != 8):
    print "Usage: mpirun python antibody_mpi.py modelname ndecoys <light-chain.fasta> <heavy-chain.fasta> antibody_path rosetta_path rosettadb_path"
    exit()
modelname = sys.argv[1].strip()
ndecoys = int(sys.argv[2])
lightchain = sys.argv[3].strip()
heavychain = sys.argv[4].strip()
antibody = sys.argv[5].strip()
rosetta = sys.argv[6].strip()
# Try to find the default, single processor Rosetta binary
platform = None
binaries = glob.glob(rosetta + "/antibody_graft.default.*")
if (len(binaries) > 0):
    platform = binaries[0].split("antibody_graft")[1]
    platform = platform[1:] # Get rid of the leading .
else:
    binaries = glob.glob(rosetta + "/antibody_graft.static.*")
    if (len(binaries) > 0):
        platform = binaries[0].split("antibody_graft")[1]
        platform = platform[1:] # Get rid of the leading
rosettadb = sys.argv[7].strip()
noderank = MPI.COMM_WORLD.Get_rank()
nodesize = MPI.COMM_WORLD.Get_size()
print "platform=" + platform
for i in range(noderank, ndecoys, nodesize):
    # Note to future developers: the --rosetta-platform argument is essential
    # If you do not include it, it will default to <rosetta_executable>.linuxgccrelease
    # This a problem, because that executable is whatever the last build type was
    # So if the last build was the MPI build (what it usually is), then that executable is the MPI one
    # If you try to run the MPI executable within a Python script that is mpi4py aware, it crashes
    if (platform):
        ABargs = "python " + antibody + "/antibody.py --light-chain " + lightchain + " --heavy-chain " + heavychain + " --antibody-database=" + antibody + "/antibody_database --blast-database=" + antibody + "/blast_database --rosetta-platform=" + platform + " --rosetta-bin=" + rosetta + " --rosetta-database=" + rosettadb + " --prefix=grafting" + str(i+1) + "/"
    else:
        ABargs = "python " + antibody + "/antibody.py --light-chain " + lightchain + " --heavy-chain " + heavychain + " --antibody-database=" + antibody + "/antibody_database --blast-database=" + antibody + "/blast_database --rosetta-bin=" + rosetta + " --rosetta-database=" + rosettadb + " --prefix=grafting" + str(i+1) + "/"
    print "Starting antibody model " + str(i+1) + "..."
    f = open("out" + str(i), "w")
    f.write(ABargs + "\n")
    f.close()
    (out, err) = subprocess.Popen(ABargs, shell=True, stdout=subprocess.PIPE).communicate()
    print ABargs
    # The Chothia numbering gives rise to duplicate residues, that BioPython cannot process
    # This gives all the residues unique numbers so they are not issues later on
    cleanPDB("grafting" + str(i+1) + "/grafted.relaxed.pdb")
    os.rename("grafting" + str(i+1) + "/grafted.relaxed.pdb", modelname + ("_%4.4i" % (i+1)) + ".pdb")
    shutil.rmtree("grafting" + str(i+1), ignore_errors=True)
MPI.Finalize()
