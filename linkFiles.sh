#!/bin/bash

SRC=$HOME/versionedScriptsAndDotfiles

mkdir -p $HOME/bin
ln -sf $SRC/generalScripts $HOME/bin/
ln -sf $SRC/biopython $HOME/bin/
ln -sf $SRC/pyscripts $HOME/bin/

exit

TARG=$HOME/.bashrc
if [  ! -h "$TARG" ];then
	if [ -e "$TARG" ];then
		mv $TARG $TARG.bak
	fi
	ln -s $HOME/$SRC/dot_profile $TARG
fi

TARG=$HOME/.vimrc
if [  ! -h "$TARG" ];then
	if [ -e "$TARG" ];then
		mv $TARG $TARG.bak
	fi
	ln -s $HOME/$SRC/dot_vimrc $TARG
fi

TARG=$HOME/.vim
if [  ! -h "$TARG" ];then
	if [ -d "$TARG" ];then
		mv $TARG $TARG.bak
	fi
	ln -s $HOME/$SRC/vim_dir $TARG
fi




