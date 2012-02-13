#!/bin/bash


FROM=${1%/}

if [ ! -d "$1" ];then
    echo $1 is not a directory!
    exit
fi

if [ ! "$#" -eq 2 ];then
    echo "WARNING - assuming same name for destination directory"
    TO=${1%/}
    TO=${TO##*/}
else
    TO=${2%/}
fi

rsync -av --include=*/ --include=*.source --include=*.py --include=*.sh --exclude=*  $FROM/ $TO/
