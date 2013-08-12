#!/bin/bash

<<TEMP
###########################
#General stuff
###########################
ssh-keygen

sudo apt-get install --assume-yes git vim
sudo apt-get install --assume-yes build-essential autoconf subversion libtool

sudo apt-get install --assume-yes matplotlib python-numpy python-scipy r-base 

git clone ssh://zwickl@phylo.bio.ku.edu/~/gitRepos/versionedScriptsAndDotfiles.git
$HOME/versionedScriptsAndDotfiles/linkFiles.sh
source $HOME/.bashrc
TEMP

<<TEMP
###########################
#NCL
###########################
mkdir -p $HOME/ncl
NROOT=$HOME/ncl/ncl-2.1
svn checkout http://svn.code.sf.net/p/ncl/code/branches/v2.1 $NROOT
mkdir -p $NROOT/builds/gccNoDynamic
cd $NROOT
./bootstrap.sh
cd builds/gccNoDynamic
env CXXFLAGS="-DNCL_CONST_FUNCS" ../../configure --prefix=`pwd`/installed --disable-shared
make && make check && make install
cd $HOME
###########################

<<TEMP
###########################
#BOINC
###########################
git clone git://boinc.berkeley.edu/boinc-v2.git boinc-v2
sudo apt-get install libcurl4-openssl-dev libfcgi-dev
cd boinc-v2
./_autosetup
./configure
exit
TEMP

<<TEMP
###########################
#GARLI
###########################
svn checkout https://garli.googlecode.com/svn/garli googleCodeRepo --username dzwickl@gmail.com
TEMP

###########################
#for beagle
###########################
sudo apt-get install nvidia-opencl-dev
sudo apt-get install default-jdk
sudo apt-get install nvidia-cuda-dev
sudo apt-get install nvidia-cuda-toolkit

svn co https://beagle-lib.googlecode.com/svn/ beagle --username=dzwickl@gmail.com


