#!/bin/sh

set -e

brew install automake autoconf libtool gpatch
ln -s $(which glibtoolize) /usr/local/bin/libtoolize
ln -s $(which gpatch) /usr/local/bin/patch

pushd /tmp
wget --no-verbose https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz
tar xzf gsl-2.7.1.tar.gz

pushd gsl-2.7.1
wget https://raw.githubusercontent.com/Homebrew/formula-patches/03cf8088210822aa2c1ab544ed58ea04c897d9c4/libtool/configure-big_sur.diff
patch -p 1 < configure-big_sur.diff
./configure --disable-dependency-tracking --prefix=/usr/local/
make
sudo make install
popd # gsl-2.7.1

git clone https://github.com/scipy/boost-headers-only.git
sudo mv boost-headers-only/boost /usr/local/include/boost

popd # /tmp
