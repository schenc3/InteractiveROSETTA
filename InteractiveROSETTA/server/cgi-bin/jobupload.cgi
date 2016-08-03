#!/usr/bin/python
import sys
import os
import commands

iRosetta_home = ".." # Alter this if iRosetta home is not the parent directory of this script
hostlist = iRosetta_home + "/hostlist"

print "Content-type:text/html\n\n"
try:
    # Get the index of the last machine to get a job, so we can send this job to the next one
    try:
        f = open("machineID", "r")
        machineID = int(f.readlines())
        f.close()
    except:
        machineID = -1
    machineID = machineID + 1
    f = open(hostlist, "r")
    machines = []
    for aline in f:
        machines.append(aline.strip())
    f.close()
    if (machineID >= len(machines)):
        machineID = 0
    worker = machines[machineID]
    filedata = []
    while (True):
        aline = sys.stdin.readline()
        if (not(aline)):
            break
        if (aline[0:19] == "Content-Disposition"):
            data = aline.split()
            for val in data:
                if (val.find("name=") >= 0 and val.find("filename=") < 0):
                    filename = val[val.find("\"")+1:val.rfind("\"")]
        elif (aline[0:2] == "--"):
            ID = aline[2:len(aline.strip())-2]
        elif (aline.find("Content") < 0):
            filedata.append(aline.strip())
    # Make sure the filename is something we'd be expecting to see from InteractiveROSETTA
    if (not(filename.endswith("input"))):
        raise Exception()
    if (filename == "testinput"):
        print "InteractiveROSETTA Upload Successful"
    else:
        # Now make sure the ID is not something potentially malicious...
        ID = ID.lower()
        for i in range(len(ID)-1, -1, -1):
            if (not(ID[i] in "abcdefghijklmnopqrstuvwxyz0123456789")):
                ID = ID[0:i] + "0" + ID[i+1:]
        f = open(iRosetta_home + "/jobfiles/uploaded.0.0.0-" + filename + "-" + ID, "w")
        for aline in filedata:
            f.write(aline + "\n")
        f.close()
        os.rename(iRosetta_home + "/jobfiles/uploaded.0.0.0-" + filename + "-" + ID, iRosetta_home + "/jobfiles/" + worker + "-" + filename + "-" + ID)
        f = open("machineID", "w")
        f.write(str(machineID))
        f.close()
        print "InteractiveROSETTA Upload Successful"
except:
    print "InteractiveROSETTA Upload Failed"
