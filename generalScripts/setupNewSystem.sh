#!/bin/bash

sudo apt-get install git vim

git clone ssh://zwickl@phylo.bio.ku.edu/~/gitRepos/versionedScriptsAndDotfiles.git

sudo apt-get install build-essential  subversion matplotlib python-numpy python-scipy r-base 

svn checkout https://garli.googlecode.com/svn/garli/ googleCodeRepo --username dzwickl@gmail.com 
svn checkout http://svn.code.sf.net/p/ncl/code/branches/v2.1 ncl-2.1

