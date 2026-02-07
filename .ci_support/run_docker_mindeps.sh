#!/bin/sh

set -xe

apt-get -y update && apt-get -y install curl gnupg

echo deb [signed-by=/usr/share/keyrings/openturns-keyring.gpg] https://openturns.github.io/apt/debian bullseye main | tee /etc/apt/sources.list.d/openturns.list
curl -fsSL https://openturns.github.io/apt/public.key | gpg --dearmor --yes --output /usr/share/keyrings/openturns-keyring.gpg

apt-get -y update && apt-get -y install git g++ cmake swig python3-dev python3-openturns libopenturns-dev libcgal-dev libnanoflann-dev

set -e
git config --global --add safe.directory /io

cd /tmp
cmake -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic -D_GLIBCXX_ASSERTIONS" -DCMAKE_INSTALL_PREFIX=$PWD/install -DSWIG_COMPILE_FLAGS="-O1" -B build /io
cd build
make install
ctest -R pyinstallcheck --output-on-failure --schedule-random ${MAKEFLAGS}
