### IMPORT WHATEVER YOU NEED
### ...
### ...
from rosetta import *
import time

### THIS FUNCTION MUST BE CALLED "runJob" AND IT ACCEPTS THE PROTOCOL KEYWORD
def runJob(protocolcode):
    ### INPUT/OUTPUT NAMES
    inputfile = protocolcode + "input"
    outputfile = protocolcode + "output"
    try:
	fin = open(inputfile, "r")
    except:
	raise Exception("ERROR: The file \"" + inputfile + "\" is missing!")
    ### PARSE THE PROTOCOL INPUT FROM THE GUI
    ### ...
    ### ...
    ### ======================================================================================
    ### THIS USES THE EXAMPLE GIVEN IN THE GUI
    ### DELETE ME LATER WHEN YOU ADD YOUR OWN CODE
    writingData = False
    for aline in fin:
	if (aline.startswith("PDBFILE")):
	    pdbfile = aline.split("\t")[1].strip()
	elif (aline.startswith("RESFILE")):
	    resfile = aline.split("\t")[1].strip()
	elif (aline.startswith("BEGIN PDB DATA")):
	    fout = open(pdbfile, "w")
	    writingData = True
	elif (aline.startswith("BEGIN RESFILE DATA")):
	    fout = open(resfile, "w")
	    writingData = True
	elif (aline.startswith("END PDB DATA") or aline.startswith("END RESFILE DATA")):
	    writingData = False
	elif (writingData):
	    fout.write(aline)
    ### END DELETION BLOCK
    ### ======================================================================================
    fin.close()
    fout = open(outputfile + "temp", "w")
    ### CODE TO RUN YOUR PROTOCOL
    ### ...
    ### ...
    ### ======================================================================================
    ### THIS EXAMPLE JUST RUNS A SIMPLE FIXBB DESIGN
    ### DELETE IT WHEN YOU WRITE YOUR OWN CODE
    init(extra_options="-ignore_unrecognized_res")
    pose = pose_from_pdb(pdbfile)
    scorefxn = create_score_function("talaris2013")
    design_pack = TaskFactory.create_packer_task(pose)
    parse_resfile(pose, design_pack, resfile)
    pack_mover = PackRotamersMover(scorefxn, design_pack)
    # You can display a progress bar in the main GUI
    # Simply create a file called "progress"
    # The first line of this file should be a fraction representing the progress (e.g. 4/5)
    # The second line is optional and contains a message that will be displayed in the dialog
    # If no message is given, it retains the previously-given message
    fout = open("progress", "w")
    fout.write("0/10\n")
    fout.write("Performing protein design...\n")
    fout.close()
    try:
	pack_mover.apply(pose)
    except Exception as e:
	raise Exception("ERROR: A problem was encountered:\n" + e.message)
    pose.dump_pdb("output.pdb")
    # Run a sleep loop just to demo the progress bar
    for i in range(1, 11):
	fout = open("progress", "w")
	fout.write(str(i) + "/10\n")
	fout.write(pose.residue(i).name3() + " was selected at position " + str(i) + "\n")
	fout.close()
	time.sleep(2)
    # Write the output file
    # It needs to have temp appended so that the file is not read by the GUI before the
    # daemon finishes writing all the output
    fout = open(outputfile + "temp", "w")
    fout.write("PDBFILE\toutput.pdb\n")
    ### END DELETION BLOCK
    ### ======================================================================================
    fout.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename(outputfile + "temp", outputfile)