# This script is used on OSX to determine if the user has installed wxPython and
# XQuartz.  It will use Tkinter to show dialogs to instruct the user what is occuring
# These pkg files will be packaged with InteractiveROSETTA so there is no need to download
# them from the Internet at any point.
from Tkinter import *
from tkMessageBox import *
import os
import sys
import commands

root = Tk()
root.geometry("1x1")

scriptdir = os.getcwd()
indx = sys.argv[0].rfind("/")
if (indx >= 0):
    tempdir = sys.argv[0][0:indx]
else:
    tempdir = ""
if (len(tempdir) > 0 and tempdir[0] == "/"):
    scriptdir = tempdir
else:
    scriptdir = scriptdir + "/" + tempdir
# Is wx installed?
try:
    import wx
except:
    answer = askyesno("wxPython Not Installed", "wxPython is not installed.  Would you like to install it now?  You will need admin privileges to continue.")
    if (answer):
        showinfo("Install wxPython", "Please follow the installation instructions")
	# res, output = commands.getstatusoutput("hdiutil attach -nobrowse " + scriptdir + "/Contents/Resources/wxPython3.0-osx-3.0.2.0-cocoa-py2.7.dmg")
 #        res, output = commands.getstatusoutput("open -W /Volumes/wxPython3.0-osx-3.0.2.0-cocoa-py2.7/wxPython3.0-osx-cocoa-py2.7.pkg")
        res, output = commands.getstatusoutput("open -W %s/Contents/Resources/wxPython3.0-osx-cocoa-py2.7.pkg"%(scriptdir))
	if (res):
	    showerror("Cannot Run InteractiveROSETTA", "InteractiveROSETTA cannot be run until wxPython is installed.")
	    print "Aborted"
            # res, output = commands.getstatusoutput("umount /Volumes/wxPython3.0-osx-3.0.2.0-cocoa-py2.7")
	    exit()
        # res, output = commands.getstatusoutput("umount /Volumes/wxPython3.0-osx-3.0.2.0-cocoa-py2.7")
    else:
	showerror("Cannot Run InteractiveROSETTA", "InteractiveROSETTA cannot be run until wxPython is installed.")
	print "Aborted"
	exit()
# Is xQuartz installed?
if (not(os.path.exists("/opt/X11"))):
    answer = askyesno("XQuartz Not Installed", "XQuartz is not installed.  Would you like to install it now?  You will need admin privileges to continue.")
    if (answer):
        showinfo("Install XQuartz", "Please follow the installation instructions")
        res, output = commands.getstatusoutput("hdiutil attach -nobrowse " + scriptdir + "/Contents/Resources/XQuartz.dmg")
        res, output = commands.getstatusoutput("open -W /Volumes/XQuartz-2.7.7/XQuartz.pkg")
	if (res):
	    showerror("Cannot Run InteractiveROSETTA", "InteractiveROSETTA cannot be run until XQuartz is installed.")
	    print "Aborted"
            res, output = commands.getstatusoutput("umount /Volumes/XQuartz-2.7.7")
	    exit()
        res, output = commands.getstatusoutput("umount /Volumes/XQuartz-2.7.7")
    else:
	showerror("Cannot Run InteractiveROSETTA", "InteractiveROSETTA cannot be run until XQuartz is installed.")
	print "Aborted"
	exit()
#this is made unecessary now that everything has been relinked
'''
#Ensure All Libraries are installed
dylibs = ['libfreetype.6.dylib','libgfortran.3.dylib','libgfortran.dylib','libGLEW.1.11.0.dylib','libopenbabel.4.dylib','libpng16.16.dylib']
for dylib in dylibs:
    res, output = commands.getstatusoutput('find /usr/local/lib -name "%s" -print'%(dylib))
    if dylib not in output:
        answer = askyesno("Library %s not found.  Install? y/n"%(dylib))
        if answer:
            showinfo('Installing library: %s'%(dylib))
            res, output = commands.getstatusoutput('cp -v %s/lib/%s /usr/local/lib/'%(scriptdir,dylib))
            if res:
                showerror("Cannot Install Library %s"%(dylib),output)
                print "Aborted"
                exit()
        else:
            showerror("Cannot Run InteractiveROSETTA","InteractiveROSETTA cannot be run until %s is installed."%(dylib))
            print "Aborted"
            exit()
'''
