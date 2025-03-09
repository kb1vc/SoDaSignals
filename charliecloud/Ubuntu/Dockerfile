FROM ubuntu:latest

RUN apt-get update -y
RUN apt-get install -y build-essential
RUN apt-get install -y libfftw3-dev
RUN apt-get install -y git cmake make g++

RUN mkdir /soda_testbuild
WORKDIR /soda_testbuild
RUN git clone https://github.com/kb1vc/SoDaSignals.git
RUN git clone https://github.com/kb1vc/SoDaUtils.git
WORKDIR /soda_testbuild/SoDaSignals
RUN mkdir build
WORKDIR /soda_testbuild/SoDaSignals/build
RUN cmake ../
RUN make
RUN ctest
RUN make install
WORKDIR /soda_testbuild/SoDaUtils
RUN mkdir build
WORKDIR /soda_testbuild/SoDaUtils/build
RUN cmake ../
RUN make
RUN ctest -j 4
RUN make install
RUN mkdir /soda_testbuild/example_test
WORKDIR /soda_testbuild/example_test
RUN cp /usr/local/share/SoDaSignals/examples/* .
RUN mkdir build
WORKDIR  /soda_testbuild/example_test/build
RUN cmake ../
RUN make
