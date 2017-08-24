import commands
import sys
import shutil

if (len(sys.argv) != 2):
    print "Usage: python killjob.py jobID"
    exit()
if (sys.argv[0].rfind("/") >= 0):
    scriptdir = sys.argv[0][0:sys.argv[0].rfind("/")]
else:
    scriptdir = "."
res, psoutput = commands.getstatusoutput("ps aux | grep " + sys.argv[1].strip() + " | grep -v grep | grep -v killjob")
# Now kill each process with this jobID
for process in psoutput.split("\n"):
    if (len(process.strip()) == 0):
        continue
    pid = process.split()[1]
    res, output = commands.getstatusoutput("kill -9 " + pid)
shutil.rmtree(scriptdir + "/results/" + sys.argv[1].strip(), ignore_errors=True)