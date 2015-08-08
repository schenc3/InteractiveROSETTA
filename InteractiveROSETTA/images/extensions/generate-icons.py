import commands
import glob

iconsizes = [8, 16, 22, 24, 32, 48, 64, 96, 256]
for PNG in glob.glob("*.png"):
    # Ignore resized icons
    if ("-" in PNG):
	continue
    for size in iconsizes:
	commandline = "convert " + PNG + " -resize " + str(size) + "x" + str(size) + " " + PNG.split(".png")[0] + "-" + str(size) + ".png"
	res, output = commands.getstatusoutput(commandline)