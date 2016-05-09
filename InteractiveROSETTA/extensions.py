import sys
import os.path
import platform

iRosetta_path = sys.argv[1]
icons = [("cfg", "Rosetta Configuration File"),
         ("cst", "Constraints File"),
         ("ensb", "Ensemble Archive"),
         ("fasta", "FASTA Sequence"),
         ("gz", "Compressed GZ Archive"),
         ("irm", "InteractiveROSETTA Module"),
         ("msdar", "Multi-State Design Archive"),
         ("msd", "Multi-State Design Save"),
         ("params", "Rosetta Parameters"),
         ("pdb", "ProteinDataBank File"),
         ("py", "Python Script"),
         ("resfile", "Rosetta Design Palette"),
         ("scan", "Point Mutant Scan")]
if (platform.system() == "Windows"):
    # We need to update the Windows registry
    # These changes probably won't take effect until the user logs back in
    import _winreg
    for icon, desc in icons:
        # Create filetype key
        regpath = "Software\\Classes\\." + icon
        try:
            _winreg.CreateKey(_winreg.HKEY_CURRENT_USER, regpath)
        except:
            pass
        reg_key = _winreg.OpenKey(_winreg.HKEY_CURRENT_USER, regpath, 0, _winreg.KEY_WRITE)
        _winreg.SetValueEx(reg_key, "", 0, _winreg.REG_SZ, icon + "_auto_file")
        _winreg.CloseKey(reg_key)
        # Create description for file type
        regpath = "Software\\Classes\\" + icon + "_auto_file"
        try:
            _winreg.CreateKey(_winreg.HKEY_CURRENT_USER, regpath)
        except:
            pass
        reg_key = _winreg.OpenKey(_winreg.HKEY_CURRENT_USER, regpath, 0, _winreg.KEY_WRITE)
        _winreg.SetValueEx(reg_key, "", 0, _winreg.REG_SZ, desc)
        _winreg.CloseKey(reg_key)
        # Create icon for filetype
        regpath = "Software\\Classes\\" + icon + "_auto_file\\DefaultIcon"
        try:
            _winreg.CreateKey(_winreg.HKEY_CURRENT_USER, regpath)
        except:
            pass
        reg_key = _winreg.OpenKey(_winreg.HKEY_CURRENT_USER, regpath, 0, _winreg.KEY_WRITE)
        _winreg.SetValueEx(reg_key, "", 0, _winreg.REG_SZ, iRosetta_path + "\\images\\extensions\\" + icon + ".ico")
        _winreg.CloseKey(reg_key)
        if (icon == "py" or icon == "gz"):
            continue
        # Default program key
        regpath = "Software\\Classes\\" + icon + "_auto_file\\Shell\\Open\\Command"
        try:
            _winreg.CreateKey(_winreg.HKEY_CURRENT_USER, regpath)
        except:
            pass
        reg_key = _winreg.OpenKey(_winreg.HKEY_CURRENT_USER, regpath, 0, _winreg.KEY_WRITE)
        _winreg.SetValueEx(reg_key, "", 0, _winreg.REG_SZ, "\"" + iRosetta_path + "\\InteractiveROSETTA.bat\" \"%1\"")
        _winreg.CloseKey(reg_key)
elif (platform.system() == "Linux"):
    import commands
    import glob
    # Get a list of all the icons
    # This script should already have root access to be able to change these things
    iconthemes = glob.glob("/usr/share/icons/*")
    for icon, desc in icons:
        # Update the mime information
        try:
            fout = open("/usr/share/mime/packages/" + icon + ".xml", "w")
            fout.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>")
            fout.write("<mime-info xmlns=\"http://www.freedesktop.org/standards/shared-mime-info\">")
            fout.write("  <mime-type type=\"text/" + icon + "\">")
            fout.write("    <!--Created automatically by update-mime-database. DO NOT EDIT!-->")
            fout.write("    <comment>" + desc + "</comment>")
            fout.write("    <icon name=\"text-" + icon + "\"/>")
            fout.write("    <glob-deleteall/>")
            fout.write("    <glob pattern=\"*." + icon + "\"/>")
            fout.write("  </mime-type>")
            fout.write("</mime-info>")
            fout.close()
            # Now try to add all the icons to the themes
            for theme in iconthemes:
                # Get a list of the file sizes
                iconsizedirs = glob.glob(theme + "/*")
                for iconsizedir in iconsizedirs:
                    if (iconsizedir[len(iconsizedir)-1] not in "0123456789" or "x" not in iconsizedir):
                        continue
                    iconsize = iconsizedir[iconsizedir.rfind("/")+1:]
                    try:
                        size = iconsize[0:iconsize.find("x")]
                        commandline = "cp \"" + iRosetta_path + "/images/extensions/" + icon + "-" + size + ".png\" " + iconsizedir + "/mimetypes/text-" + icon + ".png"
                        res, output = commands.getstatusoutput(commandline)
                        if (res):
                            raise Exception()
                        if (icon == "pdb"):
                            # Extra stuff
                            commandline = "cp \"" + iRosetta_path + "/images/extensions/" + icon + "-" + size + ".png\" " + iconsizedir + "/mimetypes/chemical-x-" + icon + ".png"
                            res, output = commands.getstatusoutput(commandline)
                            if (res):
                                raise Exception()
                        elif (icon == "py"):
                            # Extra stuff
                            commandline = "cp \"" + iRosetta_path + "/images/extensions/" + icon + "-" + size + ".png\" " + iconsizedir + "/mimetypes/text-x-python.png"
                            res, output = commands.getstatusoutput(commandline)
                            if (res):
                                raise Exception()
                    except:
                        pass
                commandline = "gtk-update-icon-cache " + theme
                res, output = commands.getstatusoutput(commandline)
                try:
                    if (res):
                        raise Exception()
                except:
                    pass
        except:
            print "Failed to craete a MIMEtype for " + icon
    commandline = "update-mime-database /usr/share/mime"
    res, output = commands.getstatusoutput(commandline)
