# SoDaSignals Simple signal processing classes

 To get to doxygen generated documentation that is more detailed 
 than the following summary, go [here](https://kb1vc.github.io/SoDaSignals/)

 SoDaSignals contains classes that provide:
 . an interface to FFTW for complex DFT and IDFT
 . a filter class for complex signals
 . a digital down converter and up converter
 
## Dependencies

SoDaSignals requires very little in terms of libraries and other
stuff.
  . FFTW3: both double and single precision versions. If your
  distribution distinguishes between "normal" kits and "devel" kits,
  install the "devel" package. 
  . C++: Any version supporting the C++11 standard

and that's it.

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
and the include in /usr/local/include/SoDaFormat/Format.hxx

It will also write doxygen output that starts at /usr/local/share/sodaformat/doc/html/index.html

### MacOS

Mac users should probably use "/opt/local" as the install prefix. This
will put the library where macports puts other libraries and such.

For reasons I can't fully comprehend, some installations (especially
for MacOS BigSur) have trouble finding the FFTW package with cmake. As
ugly as it is, you can fix this problem like this: (in the build
directory) ``` FFTW3f_DIR=/opt/local/lib/cmake/fftw3 cmake ../ ```
Replace /opt/local/lib/cmake/fftw3 with the directory containing the
file "FFTW3Config.cmake"

Additionally, the macports FFTW3 installation is broken.  It is
missing "FFTW3fLibraryDepends.cmake" If you are stuck, you can copy
the file in "cmake/macports_fix/FFTW3LibraryDepends.cmake" to the
directory mentioned above. ("/opt/local....fftw3" or whatever)

### Linux

No trouble reported yet

### Windows

I have no idea.  


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
Did that work?

If so, congratulations.  And you can use the "CMakeLists.txt" file in
the example directory as a starting point for including the library in
your own cmake based projects.

