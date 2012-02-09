#!/bin/bash

if [ ! "$#" -eq 2 ];then
    echo "enter source and destination directory names"
    exit
fi

if [ ! -d "$1" ];then
    echo $1 is not a directory!
    exit
fi

FROM=${1%/}
TO=${2%/}

rsync -av --include=*/ --include=*.source --include=*.py --include=*.sh --exclude=*  $FROM/ $TO/
