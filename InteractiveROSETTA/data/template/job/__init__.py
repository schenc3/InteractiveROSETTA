### IMPORT WHATEVER YOU NEED
### ...
### ...
from rosetta import *

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
    fin.close()
    fout = open(outputfile + "temp", "w")
    ### CODE TO RUN YOUR PROTOCOL
    ### ...
    ### ...
    fout.close()
    # So the main GUI doesn't attempt to read the file before the daemon finishes writing its contents
    os.rename(outputfile + "temp", outputfile)