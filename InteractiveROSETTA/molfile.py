# This is a script to unpack the molfile2params.tar.gz archive so we can get the residue/ligand creator to work
# I cannot package these scripts in InteractiveROSETTA because they are licensed by RosettaCommons
# Once the user downloads PyRosetta separately, this script can unpackage the files and put them where they need to be
# This needs to run with Administrator or root privileges since InteractiveROSETTA is usually installed globally

import os
import os.path
import tarfile
import time
import platform
import glob
import wx

if (platform.system() != "Windows" and os.path.exists("scripts/rosetta_py")):
    exit()
elif (platform.system() == "Windows" and os.path.exists("scripts\\rosetta_py")):
    exit()

app = wx.App()
busyDlg = wx.BusyInfo("Searching for PyRosetta Installation...")
print "Importing molfile2params, please be patient..."
cfgfile = os.path.expanduser("~") + "/InteractiveROSETTA/seqwindow.cfg"
try:
    f = open(cfgfile.strip(), "r")
    rosettadir = "Not Found"
    rosettadb = "Not Found"
    for aline in f:
        if ("[ROSETTAPATH]" in aline):
            rosettapath = aline.split("\t")[1].strip()
        if ("[ROSETTADB]" in aline):
            rosettadb = aline.split("\t")[1].strip()
    f.close()
    if (rosettapath == "Not Found"):
        raise Exception
    if (platform.system() == "Windows"):
        molfiletgz = rosettapath + "\\toolbox\\molfile2params.tar.gz"
    else:
        molfiletgz = rosettapath + "/toolbox/molfile2params.tar.gz"
    if (not(os.isfile(molfiletgz))):
        raise Exception
except:
    # The error may have been the Rosetta import, which means the file needs to be closed
    try:
        f.close()
    except:
        pass
    # Okay, we still didn't get it, so let's traverse the filesystem looking for it...
    foundIt = False
    if (platform.system() == "Windows"):
        roots = ["C:\\"]
    elif (platform.system() == "Darwin"):
        roots = ["/Users", "/Applications", "/"]
    else:
        roots = ["/"]
    for root in roots:
        for dpath, dnames, fnames in os.walk(root):
            try:
                if (platform.system() == "Windows"):
                    try:
                        indx = fnames.index("rosetta.pyd") # 64bit
                    except:
                        indx = fnames.index("rosetta.dll") # 32bit
                else:
                    indx = dnames.index("rosetta")
                    files = glob.glob(dpath + "/rosetta/*libmini*")
                    if (len(files) == 0):
                        raise Exception
            except:
                continue
            foundIt = True
            rosettapath = dpath
            for dname in dnames:
                if ("database" in dname):
                    rosettadb = dpath + "/" + dname
                    break
            break
        if (foundIt):
            break
    if (not(foundIt)):
        print "Cannot find PyRosetta on your computer.  Did you install it?"
        time.sleep(10)
    else:
        data = []
        f = open(cfgfile, "r")
        for aline in f:
            if (not("[ROSETTAPATH]" in aline) and not("[ROSETTADB]") in aline):
                data.append(aline.strip())
        f.close()
        f = open(cfgfile, "w")
        for aline in data:
            f.write(aline + "\n")
        f.write("[ROSETTAPATH]\t" + rosettapath.strip() + "\n")
        f.write("[ROSETTADB]\t" + rosettadb.strip() + "\n")
        f.close()
busyDlg = None
if (platform.system() == "Windows"):
    molfiletgz = rosettapath + "\\toolbox\\molfile2params.tar.gz"
else:
    molfiletgz = rosettapath + "/toolbox/molfile2params.tar.gz"
os.chdir("scripts")
tar = tarfile.open(molfiletgz, "r:gz")
for item in tar:
    tar.extract(item)
tar.close()
for dpath, dnames, fnames in os.walk("molfile2params"):
    if (platform.system() == "Windows"):
        newpath = dpath[dpath.find("\\")+1:]
    else:
        newpath = dpath[dpath.find("/")+1:]
    if (newpath == dpath):
        newpath = "."
    for dname in dnames:
        if (not(os.path.exists(newpath + "/" + dname))):
            os.mkdir(newpath + "/" + dname)
    for fname in fnames:
        os.rename(dpath + "/" + fname, newpath + "/" + fname)
app.MainLoop()