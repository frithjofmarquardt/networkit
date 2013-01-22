import os
import fnmatch

# SOURCE
source = []

# walk source directory and find .cpp and .h
for (dirpath, dirnames, filenames) in os.walk("src"):
    for name in fnmatch.filter(filenames, "*.h"):
        source.append(os.path.join(dirpath, name))
    for name in fnmatch.filter(filenames, "*.cpp"):
        source.append(os.path.join(dirpath, name))


# exclude files matching following patterns
xpatterns = ["*ParametrizedGTest.*"]
excluded = []

for pattern in xpatterns:
	for name in fnmatch.filter(source, pattern):
		excluded.append(name)

source = [name for name in source if name not in excluded]


# ENVIRONMENT

## environment: macbook
macbook = Environment()
### include
macbook.Append(CPPPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
                          "~/workspace/gtest/include", \
                          "/usr/local/Cellar/log4cxx/0.10.0/include", \
                          "~/workspace/STINGER/include"])
macbook.Append(CPATH = ["~/workspace/STINGER/include"])

print(macbook["CPPPATH"])
print(macbook["CPATH"])

### link
macbook.Append(LIBS = ["STINGER", "gtest", "log4cxx"])
macbook.Append(LIBPATH = ["~/workspace/STINGER/OpenMP Debug",\
                           "~/workspace/gtest/lib", \
                            "/usr/local/Cellar/log4cxx/0.10.0/lib"])

### compiler & flags
macbook["CXX"] = "g++-4.7"
macbook.Append(CCFLAGS = "-O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11")
# TODO: extract environment-independent flags


## environment: compute11

compute11 = Environment()
compute11.Append(LIBS = ["STINGER", "gtest", "log4cxx"])


# TODO: for gcc-4.6 env.Append(CCFLAGS = "-O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11")


## select environment
# custom command line options
AddOption("--machine",
          dest="machine",
          type="string",
          nargs=1,
          action="store",
          help="specify the machine (environment) on which to build")


environments = {"macbook" : macbook, "compute11" : compute11}

try:
    env = environments[GetOption("machine")]
except:
    print("ERROR: In order to build call scons with --machine=<MACHINE> where <MACHINE> is one of: %s" % environments.keys())
    exit()


# TARGET
env.Program("EnsembleClustering-DPar", source)