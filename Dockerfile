FROM ubuntu:25.10

RUN apt-get update && apt-get install -y  \
     python3-pip\ 
     libvtk9-dev\
     cmake\
     doxygen\
     clang-20 \
     clang-tidy-20\
     clang-format-20\
     build-essential\
     wget

RUN pip install cmake-format --break-system-packages

RUN apt-get install -y gcc-15 g++-15

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-15 60 \
              --slave /usr/bin/g++ g++ /usr/bin/g++-15 && \
     update-alternatives --install /usr/bin/clang clang /usr/bin/clang-20 1000 &&\
     update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-20 1000 && \
     update-alternatives --install /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-20 1000 && \
     update-alternatives --install /usr/bin/clang-format clang-format /usr/bin/clang-format-20 1000

RUN update-alternatives --display clang && clang++ --version && clang-tidy --version && \
      g++ --version && cmake --version
