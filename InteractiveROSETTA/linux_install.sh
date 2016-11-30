#!/bin/bash

if [ -t 1 ]; then
    echo "Installing InteractiveROSETTA"
    echo ""
else
    notify-send "Please run this installer in a terminal"
    exit
fi

if [ $EUID == 0 ]; then
    echo "Do not run this script as root.  You will be prompted for the root password when necessary."
    exit
fi

# Figure out what Linux distro we are on, so we know what the package manager is
if [ -f /etc/debian_version ]; then
    PMGR="sudo apt-get -q -y install "
    ROSETTA_VER="Ubuntu"
elif [ -f /etc/redhat-release ]; then
    PMGR="sudo yum -q -y install "
    ROSETTA_VER="Scientific Linux/RedHat"
else
    PMGR="????"
    echo "Unfortunately PyRosetta is currently only supported on RedHat or Debian-based distros."
    echo "Aborting installation"
    exit
fi

REALOLDDIR="`pwd`"
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd "$SCRIPTDIR"
find ../InteractiveROSETTA -type d -exec chmod 755 {} \;
find ../InteractiveROSETTA -type f -exec chmod 755 {} \;

# Have the user unpackage PyRosetta before doing anything else
# echo "Searching for PyRosetta installation..."
# while [ 1 ]; do
#     SEARCHDIR="/"
#     read -p "Enter a directory to search (default: "$SEARCHDIR", type \"q\" to quit): " input
#     if [ $input ]; then
# 	SEARCHDIR=$input
#     fi
#     if [ $SEARCHDIR == "q" ]; then
# 	echo " "
# 	echo "Please download the PyRosetta package and unpack it to continue"
# 	echo "The PyRosetta version you need is "$ROSETTA_VER
# 	exit
#     fi
#     NEWROSETTA=`find $SEARCHDIR -name SetPyRosettaEnvironment* 2> /dev/null | head -n 1 | awk -F "/SetPyRosetta" '{print $1}' `
#     if [ $NEWROSETTA ]; then
# 	echo "Looking for PyRosetta in "$NEWROSETTA
# 	# This is obnoxious...there is a colon in the Ubuntu PyRosetta directory name
# 	# and it messes up : delimited paths, so create an appropriate symlink
# 	#if [[ $ROSETTA_VER == "Ubuntu" ]]; then
# 	    #OLDROSETTA=$NEWROSETTA
# 	    #NEWROSETTA=`echo $OLDROSETTA | awk -F "/SetPyRosetta" '{print $1}' | awk -F ":" '{print $1}'`
# 	    #sudo unlink $NEWROSETTA
# 	    #sudo ln -s $OLDROSETTA $NEWROSETTA
# 	    # Move the so files to rosetta because it cannot find them in PyRosetta root
# 	    #sudo mv $NEWROSETTA/*.so* $NEWROSETTA/rosetta
# 	#fi
# 	echo "" >> ~/.bashrc
# 	echo "source $NEWROSETTA/SetPyRosettaEnvironment.sh" >> ~/.bashrc
# 	OLDDIR=`pwd`
# 	source ~/.bashrc > /dev/null 2> /dev/null
# 	cd "$OLDDIR"
# 	break
#     else
# 	echo "A PyRosetta installation was not detected on that path!"
# 	echo "Please try again"
# 	echo " "
#     fi
# done

# Is Python installed? (Probably, but maybe not)
if hash python 2> /dev/null; then
    echo "Python installation detected"
else
    echo "Installing Python 2.7"
    $PMGR"python" > /dev/null 2> /dev/null
fi

PYTHON_MAJOR=`python -c "import sys; print sys.version_info[0]"` > /dev/null 2> /dev/null
PYTHON_MINOR=`python -c "import sys; print sys.version_info[1]"` > /dev/null 2> /dev/null
if [[ $PYTHON_MAJOR != 2 ]]; then
    echo "Python 2 is not detected.  Installing it..."
    $PMGR"python" > /dev/null 2> /dev/null
fi
if [[ $PYTHON_MINOR != 7 ]]; then
    echo "Python 2.7 is not detected.  InteractiveROSETTA has not been tested on versions of Python earlier than 2.7"
    echo "Continuing to install, but be warned that you may encounter some unexpected behavior"
fi

$PMGR"xterm" > /dev/null 2> /dev/null
echo "Installing wxPython..."
if [[ $ROSETTA_VER == "Ubuntu" ]]; then
    $PMGR"python-wxgtk2.8" > /dev/null 2> /dev/null
else
    $PMGR"wxPython" > /dev/null 2> /dev/null
fi
echo "Installing setuptools..."
$PMGR"python-setuptools" > /dev/null 2> /dev/null
$PMGR"python-dev" > /dev/null 2> /dev/null
echo "Installing PyMOL..."
#if [ -f /etc/redhat-release ]; then
    # Version 1.7.6 has some issues...use 1.7.2 instead
#    if $PMGR"`yum --showduplicates search pymol | grep pymol-1.7.2 | awk '{print $1}'`" > /dev/null 2> /dev/null; then
#	echo "PyMOL 1.7.2 installed"
#    else
#	$PMGR"pymol" > /dev/null 2> /dev/null
#	echo "PyMOL version 1.7.2 not found, installing the default instead"
#	echo "If selections are not being synchronized properly, please install v1.7.2"
#    fi
#else
$PMGR"pymol" > /dev/null 2> /dev/null
#fi
echo "Installing BioPython..."
$PMGR"python-biopython" > /dev/null 2> /dev/null
if [[ $ROSETTA_VER == "Ubuntu" ]]; then
    $PMGR"python-dev" > /dev/null 2> /dev/null
else
    $PMGR"python-devel" > /dev/null 2> /dev/null
fi
echo "Installing psutil..."
sudo easy_install -q psutil > /dev/null 2> /dev/null
echo "Installing poster..."
sudo easy_install -q poster > /dev/null 2> /dev/null
echo "Installing requests..."
sudo easy_install -q requests > /dev/null 2> /dev/null
echo "Installing pyperclip..."
sudo easy_install -q pyperclip > /dev/null 2> /dev/null
echo "Downloading OpenBabel..."
if [[ $ROSETTA_VER == "Ubuntu" ]]; then
    $PMGR"python-openbabel" > /dev/null 2> /dev/null
else
    $PMGR"openbabel" > /dev/null 2> /dev/null
    $PMGR"python-openbabel" > /dev/null 2> /dev/null
fi

# Unpack molfile2params if it has not been done yet (i.e. this is the first time this script has been run)
if [ ! -e "$SCRIPTDIR/scripts/rosetta_py" ]; then
    cp $NEWROSETTA/toolbox/molfile2params.tar.gz $SCRIPTDIR/scripts
    tar zxvf $SCRIPTDIR/scripts/molfile2params.tar.gz > /dev/null 2> /dev/null # This throws an error even though it appears to do everything right
    rm $SCRIPTDIR/scripts/molfile2params.tar.gz
    mv molfile2params/* $SCRIPTDIR/scripts
    rm -rf molfile2params
fi

# Now run a simple PyRosetta script as root so the rotamer binary files get written
echo "Generating PyRosetta binary files..."
OLDDIR=`pwd`
cd "$NEWROSETTA"
sudo cp "$OLDDIR/data/bigPDB.pdb" .
sudo python -c "from rosetta import *; init(extra_options='-ignore_unrecognized_res'); pose=pose_from_pdb('bigPDB.pdb'); scorefxn=create_score_function('talaris2013'); scorefxn(pose);" > /dev/null 2> /dev/null
sudo rm bigPDB.pdb
cd "$REALOLDDIR"

# Put the files in a global directory
echo "Installing InteractiveROSETTA to /usr/local"
sudo cp -R "$SCRIPTDIR" /usr/local
echo "InteractiveROSETTA was installed to /usr/local/InteractiveROSETTA!"
echo "Run it by executing: python /usr/local/InteractiveROSETTA/InteractiveROSETTA.py"
cd /usr/local/InteractiveROSETTA

# Update the PyMOL splash screen
SPLASHDIR=`python -c "import pymol; print pymol.__file__[0:pymol.__file__.rfind(\"/\")]"` > /dev/null 2> /dev/null
sudo mv $SPLASHDIR/pymol_path/data/pymol/splash.png $SPLASHDIR/pymol_path/data/pymol/original_splash.png
sudo cp images/pymol_splash.png $SPLASHDIR/pymol_path/data/pymol/splash.png
sudo mv $SPLASHDIR/data/pymol/splash.png $SPLASHDIR/data/pymol/original_splash.png
sudo cp images/pymol_splash.png $SPLASHDIR/data/pymol/splash.png

cd "$REALOLDDIR"

# Create a desktop shortcut
echo "[Desktop Entry]" > ~/Desktop/InteractiveROSETTA.desktop
echo "Comment=Interactive modeling of proteins using PyRosetta" >> ~/Desktop/InteractiveROSETTA.desktop
echo "Terminal=true" >> ~/Desktop/InteractiveROSETTA.desktop
echo "Name=InteractiveROSETTA" >> ~/Desktop/InteractiveROSETTA.desktop
echo "Exec=python /usr/local/InteractiveROSETTA/InteractiveROSETTA.py" >> ~/Desktop/InteractiveROSETTA.desktop
echo "Type=Application" >> ~/Desktop/InteractiveROSETTA.desktop
echo "Icon=/usr/local/InteractiveROSETTA/images/icon.png" >> ~/Desktop/InteractiveROSETTA.desktop
echo "Categories=Science;Education;" >> ~/Desktop/InteractiveROSETTA.desktop

# Try to update the MIME database to have fancy icons for the file extensions
echo "Updating MIME database..."
cd "$SCRIPTDIR"
sudo python extensions.py "$SCRIPTDIR" > /dev/null 2> /dev/null
cd "$REALOLDDIR"
