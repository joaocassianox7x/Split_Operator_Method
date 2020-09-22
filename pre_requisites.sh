#!/bin/bash

mkdir figs
mkdir out

sudo apt-get update -y

sudo apt-get install ffmpeg

#c++ libraries
sudo apt-get install -y libomp-dev
sudo apt-get install -y libblas-dev
sudo apt-get install -y libfftw3-dev
sudo apt-get install -y liblapack-dev
sudo apt-get install -ylibopenblas-dev
sudo apt-get install -y liblapacke-dev

#pytho and libraries
sudo apt-get install python3
sudo apt-get install python3-pip
pip3 install seaborn
pip3 install numpy
pip3 install scipy
pip3 install joblib
pip3 install matplotlib


#armadillo library for c++
wget http://sourceforge.net/projects/arma/files/armadillo-9.880.1.tar.xz

tar -xvf armadillo-9.880.1.tar.gz

cd armadillo-9.880.1

./configure

make
sudo make install



