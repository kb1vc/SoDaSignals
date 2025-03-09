# SoDaSignals Simple signal processing classes

 To get to doxygen generated documentation that is more detailed 
 than the following summary, build the library, install it, and
 go to {whateveryourinstallationrootis}/share/SoDaUtils/doc/html/index.html

 SoDaSignals contains classes that provide:
 . an interface to FFTW for complex DFT and IDFT
 . a filter class for complex signals
 . a digital down converter and up converter
 
## Dependencies

SoDaSignals requires very little in terms of libraries and other
stuff.
  . git -- optional
  . cmake
  . FFTW3: both double and single precision versions. If your
  distribution distinguishes between "normal" kits and "devel" kits,
  install the "devel" package. 
  . C++: Any version supporting the C++11 standard
  
and that's it.

For Fedora 
```
dnf install gcc-c++ n
dnf install git cmake
dnf install fftw-devel
```

Optionally, to get nice documentation in html form
```
dnf install doxygen
```

### For Unbuntu

```
apt install g++
apt install git
apt install cmake
apt install fftw3-dev
```

And the optional doxygen kit: 
```
apt install doxygen
```

SoDaSignals also requires the SoDaUtils library. It will attempt to
build the library itself if it isn't already installed. 

## Installing

Just like any other CMake project.  For instance, to install the
package in /usr/local ... From this directory

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
```

This will install the libraries in /usr/local/lib or lib64 as appropriate
and the includes in /usr/local/include/SoDa

It will also write doxygen output that starts at /usr/local/share/SoDaSignals/

### Installing Without root
It is possible to install the library in a private directory (without needing root) like this: 

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/my_tools/ ../
make
sudo make install
```


This should work just fine, but if you do, then any build that *uses*
SoDaSignals needs to add this to its cmake

```
cmake -DCMAKE_PREFIX_PATH=${HOME}/my_tools
```

That will tell cmake to look in your directory for the relevant cmake
files that describe where to find the libraries and headers.

### MacOS

Mac users are on their own now.... 

### Linux

No trouble reported yet

### Windows

No port exists. No port is likely. 


## Testing and Using it all

Take a look at the CMakeLists.txt file and FilterExample.cxx in the
example directory.  If the installation has gone right, then you
should be able to do this from this directory.

```
cd example
mkdir build
cd build
cmake ../
make
./FilterExample
```

If you've installed it locally in a directory like ~/foo , you'll need to do this:
```
cd example
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=${HOME}/my_tools ../
make
./FilterExample
```


Did that work?

If so, congratulations.  And you can use the "CMakeLists.txt" file in
the example directory as a starting point for including the library in
your own cmake based projects.

