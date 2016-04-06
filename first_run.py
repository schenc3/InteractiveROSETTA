#!/usr/bin/python

#Imports
import os
import os.path
import shutil
import sys
import __main__
import glob
import commands


def goToSandbox(extra=""):
    '''Easily gets us back to the sandbox location'''
    print "goToSandbox"
    homedir = os.path.expanduser("~")
    if (len(extra) > 0):
        extra = "/" + extra
    os.chdir(homedir + "/.InteractiveROSETTA" + extra)

def make_cfg(options = {'first_run':'1','lib_dir':'0'}):
    '''Generates a new osx.config file and populates it with the default options'''
    print "make_cfg:: options: ",options
    goToSandbox()
    cfg = open("osx.config",'w+')
    for option in options:
        cfg.write("%s %s\n"%(option,options[option]))
        print "%s %s"%(option,options[option])
    cfg.close()

def first_run(options,appdir=''):
    '''Does all the setup for the first run of InteractiveROSETTA'''
    print "first_run"
    libdir = "%s/lib"%(appdir)
    goToSandbox()
    #Relink all libraries to depend libraries stored in libdir
    relink_libraries(appdir,libdir)
    options['libdir'] = libdir
    options['first_run'] = '0'
    #save changes to configuration
    make_cfg(options)

def parse_otool(otool_out,lib,libdir):
    '''Given output of find otool -L and a library, determine if any libraries
    need to be relinked.  Outputs a list of tuples in the form (library,oldlib)'''
    print "parse_otool"
    otool_out = otool_out.split('\n')
    print otool_out
    needs_relink = []
    library = ''
    for line in otool_out:
        print line
        if line[0] == '\t':
            if lib in line:
                oldlib = line.split("(")[0].strip()
                print oldlib
                if oldlib != "%s/%s"%(libdir,lib):
                    needs_relink.append((library,oldlib))
                    print (library,oldlib)
        else:
            library = line.strip()


def relink_libraries(appdir,libdir):
    '''for each library in libdir, make sure that everything that points to that
    library, points to the one located in libdir'''
    print "relink_libraries"
    libs = commands.getstatusoutput("ls %s"%(libdir))[1].split()
    for lib in libs:
        #for .so files
        cmd = 'find %s -name "*.so" -exec otool -L {} ";" -or -name "*.dylib" -exec otool -L {} ";"'%(appdir)
        print cmd
        status,output = commands.getstatusoutput(cmd)
        print status,output
        # sys.exit()
        libs_to_relink = []
        if status != 0:
            raise RuntimeError("%s: %s"%(status,output))
        else:
            libs_to_relink = parse_otool(output,lib,libdir)
            print "libs_to_relink:",libs_to_relink
        if libs_to_relink != None:
            for library,oldlib in libs_to_relink:
                cmd = "install_name_tool -change %s %s/%s %s"%(oldlib,libdir,lib,library)
                print cmd
                status,output = commands.getstatusoutput(cmd)
                if status != 0:
                    raise RuntimeError("%s: %s"%(status,output))
def main():
    #save app top directory
    #Should be something along /Applications/InteractiveROSETTA.app
    #This script will live in the top directory of InteractiveROSETTA.app
    #Doing it this way accounts for any changes made to the name of the .app package
    appdir = os.path.dirname(os.path.realpath(__file__))
    #check iR directory for osx.config file
    goToSandbox()
    try:
        cfg = open("osx.config",'r')
    except IOError:
        #If osx.config doesn't exist, create it and set it to defaults.
        #Open the file now and continue as usual
        make_cfg()
        cfg = open("osx.config",'r')
    options = {}
    for line in cfg:
        line = line.split()
        options[line[0]] = ' '.join(line[1:]).strip()
    cfg.close()
    print "options:",options
    try:
        if options['first_run'] == '1':
            first_run(options,appdir)
    except KeyError:
            options["first_run"] = '1'
            first_run(options,appdir)

if __name__ == "__main__":
    main()
