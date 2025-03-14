FROM ubuntu_base

RUN rm -rf /soda_testbuild 
RUN mkdir /soda_testbuild
WORKDIR /soda_testbuild
RUN git clone https://github.com/kb1vc/SoDaSignals.git
WORKDIR /soda_testbuild/SoDaSignals
RUN git checkout ccloud_test
RUN git status
RUN mkdir build
WORKDIR /soda_testbuild/SoDaSignals/build
RUN cmake ../
RUN make
RUN ctest
RUN make install
RUN mkdir /soda_testbuild/example_test
WORKDIR /soda_testbuild/example_test
RUN cp /usr/local/share/SoDaSignals/examples/* .
RUN mkdir build
WORKDIR  /soda_testbuild/example_test/build
RUN cmake ../
RUN make
