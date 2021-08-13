# msfitslib_core
Library with classes to manipulate FITS files (including I/O). It is the same as msfitslib, but does not include general functionality classes from mscommonlib. This library needs to be linked separately

# required :
   sudo apt-get install libnova-dev libnova fftw2 fftw-dev libnova-0.16-0 libnova-dev fftw3-dev libhdf5-dev


CERN root :   
   https://root.cern.ch/building-root
   https://root.cern.ch/downloading-root

   Can be removed by editing CMakeList.txt and comment out lines :
   
   list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
   
   and
   
   find_package(ROOT REQUIRED COMPONENTS RIO Net)


# installation:

mkdir build

cd build

cmake ..

make

sudo make install

# possible problems :

   May need addition of -ldl to all target_link_libraries (if does not link) - requires to repeat steps: cmake ../;make 
