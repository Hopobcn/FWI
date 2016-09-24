#!/usr/bin/env sh


# install newer CMAKE
CLANG_TAR=cmake-3.4.3-Linux-x86_64

wget http://www.cmake.org/files/v3.4/$CLANG_TAR.tar.gz --no-check-certificate
mkdir -p $HOME/dep/cmake
tar -C $HOME/dep/cmake -xzf $CLANG_TAR.tar.gz

echo $PATH
export PATH=$HOME/dep/cmake/$CLANG_TAR/bin:$PATH
echo $PATH

