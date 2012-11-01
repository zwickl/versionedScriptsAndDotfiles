#!/bin/bash

function findcopy(){
    FILE=${1##*/}
    rsync -av --include="*/" --include="$FILE" --exclude="*" $2/ $3/
}

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

#find and copy symlinks too
#find $FROM -type l -exec /bin/cp -rp --parents '{}' $TO \;
#find $FROM -type l -exec /usr/bin/rsync -av --include="*/" --include='{}' --exclude="*" $FROM/ $TO/ \;
#find $FROM -type l | while read i ; do findcopy $i $FROM $TO ;done
