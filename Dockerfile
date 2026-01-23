FROM ubuntu:25.10

RUN apt-get update && apt-get install -y  \
     python3-pip\ 
     libvtk9-dev\
     cmake\
     doxygen\
     clang-19 \
     clang-tidy-19\
     clang-format-19\
     build-essential\
     git\
     wget

RUN pip install cmake-format --break-system-packages

RUN apt-get install -y gcc-14 g++-14

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 60 \
              --slave /usr/bin/g++ g++ /usr/bin/g++-14 && \
     update-alternatives --install /usr/bin/clang clang /usr/bin/clang-19 1000 &&\
     update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-19 1000 && \
     update-alternatives --install /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-19 1000 && \
     update-alternatives --install /usr/bin/clang-format clang-format /usr/bin/clang-format-19 1000

RUN update-alternatives --display clang && clang++ --version && clang-tidy --version && \
      g++ --version && cmake --version
