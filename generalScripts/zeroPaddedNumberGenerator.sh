#!/bin/bash

function zeroPadNumber(){
    if [ $# -gt "1" ];then
        DIGITS=$2
    else
        DIGITS=2
    fi

    RESULT=$1
    while [ ${#RESULT} -lt "$DIGITS" ]
    do
        RESULT=0$RESULT
    done
    echo $RESULT
}

SERIES=`for ((VAL=$1;VAL<=$2;VAL++))
do
    zeroPadNumber $VAL $3
done
`
echo $SERIES


