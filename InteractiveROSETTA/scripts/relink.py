from commands import getstatusoutput as run
import glob

#MUST BE RUN from the scripts directory

def link_lib(lib_file,libs):
    link_libs = [line.split()[0].strip() for line in run('otool -L %s'%(lib_file))[1].split('\n')]
    link_libs.pop(0)
    # print link_libs

    lib_libs = {}
    for lib in libs:
        for llib in link_libs:
            if lib in llib:
                lib_libs[lib] = llib
                break
    # print lib_libs

    for lib in sorted(lib_libs.keys()):
        print lib
        cmd = 'install_name_tool -change %s `pwd`/../lib/%s %s'%(lib_libs[lib], lib,lib_file)
        print cmd
        s,o = run(cmd)
        if s: print o
    s,o = run('otool -L %s'%(lib_file))
    if s: print o

#relink libraries in lib directory

libs = ['libGLEW.1.11.0.dylib','libfreetype.6.dylib','libgcc_s.1.dylib','libgfortran.3.dylib','libgfortran.dylib','libopenbabel.4.dylib','libpng16.16.dylib','libquadmath.0.dylib']
link_lib("hmmstr4py_darwin.so",libs)
other_libs = glob.glob("../lib/*.dylib")
for lib in other_libs: link_lib(lib,libs)
