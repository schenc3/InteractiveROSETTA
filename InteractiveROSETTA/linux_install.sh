#!/bin/bash

if [ $EUID == 0 ]; then
    echo "Do not run this script as root.  You will be prompted for the root password when necessary."
    exit
fi

SEARCHDIR="/"

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

# Have the user unpackage PyRosetta before doing anything else
echo "Searching for PyRosetta installation..."
read -p "Enter a directory to search (default: "$SEARCHDIR"): " input
if [ $input != "" ]; then
    SEARCHDIR=$input
fi
NEWROSETTA=`find $SEARCHDIR -name SetPyRosettaEnvironment* 2> /dev/null | head -n 1 | awk -F "/SetPyRosetta" '{print $1}' `
if [ $NEWROSETTA != "" ]; then
    # This is obnoxious...there is a colon in the Ubuntu PyRosetta directory name
    # and it messes up : delimited paths, so create an appropriate symlink
    if [[ $ROSETTA_VER == "Ubuntu" ]]; then
	OLDROSETTA=$NEWROSETTA
	NEWROSETTA=`echo $OLDROSETTA | awk -F "/SetPyRosetta" '{print $1}' | awk -F ":" '{print $1}'`
	sudo unlink $NEWROSETTA
	sudo ln -s $OLDROSETTA $NEWROSETTA
	# Move the so files to rosetta because it cannot find them in PyRosetta root
	sudo mv $NEWROSETTA/*.so* $NEWROSETTA/rosetta
    fi
    echo "" >> ~/.bashrc
    echo "source $NEWROSETTA/SetPyRosettaEnvironment.sh" >> ~/.bashrc
    OLDDIR=`pwd`
    source ~/.bashrc
    cd "$OLDDIR"
else
    echo "A PyRosetta installation was not detected."
    echo "Please download the package and unpack it."
    echo "The PyRosetta version you need is "$ROSETTA_VER
    exit
fi

# Is Python installed? (Probably, but maybe not)
if hash python 2> /dev/null; then
    echo "Python installation detected"
else
    echo "Installing Python 2.7"
    $PMGR"python"
fi

PYTHON_MAJOR=`python -c "import sys; print sys.version_info[0]";`
PYTHON_MINOR=`python -c "import sys; print sys.version_info[1]"`
if (( $PYTHON_MAJOR != 2 )); then
    echo "Python 2 is not detected.  Installing it..."
    $PMGR"python"
fi
if (( $PYTHON_MINOR != 7 )); then
    echo "Python 2.7 is not detected.  InteractiveROSETTA has not been tested on versions of Python earlier than 2.7"
    echo "Continuing to install, but be warned that you may encounter some unexpected behavior"
fi

echo "Installing wxPython..."
if [[ $ROSETTA_VER == "Ubuntu" ]]; then
    $PMGR"python-wxgtk2.8"
else
    $PMGR"wxPython"
fi
echo "Installing setuptools..."
$PMGR"python-setuptools"
$PMGR"python-dev"
echo "Installing PyMOL..."
$PMGR"pymol"
echo "Installing BioPython..."
$PMGR"python-biopython"
if [[ $ROSETTA_VER == "Ubuntu" ]]; then
	$PMGR"python-dev"
else
	$PMGR"python-devel"
fi
echo "Installing psutil..."
sudo easy_install -q psutil
echo "Installing poster..."
sudo easy_install -q poster
echo "Installing requests..."
sudo easy_install -q requests
echo "Installing OpenBabel..."
if [[ $ROSETTA_VER == "Ubuntu" ]]; then
	$PMGR"python-openbabel"
else
	$PMGR"openbabel"
	$PMGR"python-openbabel"
fi

# Unpack molfile2params if it has not been done yet (i.e. this is the first time this script has been run)
if [ ! -e "$SCRIPTDIR/scripts/rosetta_py" ]; then
    cp $ROSETTADIR/toolbox/molfile2params.tar.gz $SCRIPTDIR/scripts
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
sudo python -c "from rosetta import *; init(extra_options='-ignore_unrecognized_res'); pose=pose_from_pdb('bigPDB.pdb'); scorefxn=create_score_function('talaris2013'); scorefxn(pose);"
sudo rm bigPDB.pdb
cd "$REALOLDDIR"

# Put the files in a global directory
echo "Installing InteractiveROSETTA to /usr/local"
sudo cp -R "$SCRIPTDIR" /usr/local
echo "InteractiveROSETTA was installed to /usr/local/InteractiveROSETTA!"
echo "Run it by executing: python /usr/local/InteractiveROSETTA/InteractiveROSETTA.py"
cd /usr/local/InteractiveROSETTA

# Get the molfile stuff
sudo python molfile.py
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
