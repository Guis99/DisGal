import platform
from os.path import join

def getExecutableName(targetName, buildDir="build"):
    osName = platform.system()
    if osName == "Darwin" or osName.startswith("Linux"):
        toRun = join(".", buildDir, targetName)
    elif osName == "Windows":
        toRun = join(".", buildDir, targetName+".exe")

    return toRun