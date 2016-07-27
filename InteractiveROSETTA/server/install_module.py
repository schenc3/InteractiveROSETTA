#!/usr/bin/env python
import sys
import os
import shutil
import glob
from zipfile import ZipFile

if (len(sys.argv) != 2):
    print "Usage: python install_module.py <package.irm>"
    print "       <package.irm>: The module that will be installed"
    exit()
# Figure out where this script is running from
if ("/" in sys.argv[0]):
    home = sys.argv[0][0:sys.argv[0].rfind("/")]
else:
    home = "."
packagecode = sys.argv[1].split(".irm")[0]
packagecode = packagecode[packagecode.rfind("/")+1:]
if (os.path.exists(home + "/modules/" + packagecode)):
    while (True):
        answer = raw_input("This module already is installed.  Overwrite it? (Y/N): ").upper().strip()
        if (answer in ["Y", "N"]):
            break
    if (answer == "N"):
        print "Aborting..."
        exit()
    # Delete it
    shutil.rmtree(home + "/modules/" + packagecode, ignore_errors=True)
if (not(os.path.exists(home + "/modules"))):
    os.mkdir(home + "/modules")
os.mkdir(home + "/modules/" + packagecode)
# Unpack the irm package to the ~/InteractiveROSETTA/modules directory
fin = open(sys.argv[1].strip(), 'rb')
z = ZipFile(fin)
for name in z.namelist():
    outpath = home + "/modules/" + packagecode
    z.extract(name, outpath)
fin.close()
# Remove everything that is not the server
os.chdir(home + "/modules/" + packagecode)
files = glob.glob("*")
for filename in files:
    if (filename.strip() != "server"):
        try:
            os.remove(filename)
        except:
            try:
                shutil.rmtree(filename, ignore_errors=True)
            except:
                pass
# Now move everything up out of the server folder
for filename in glob.glob("server/*"):
    os.rename(filename, filename[filename.rfind("/")+1:])
shutil.rmtree("server", ignore_errors=True)