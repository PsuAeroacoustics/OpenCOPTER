#!/usr/bin/env bash

Xvfb :1 &> /dev/null &
xvfb_pid=$!

sleep 2

if [ $# = 1 ]; then
    DISPLAY=:1 xfoil < $1 > $2
else
    DISPLAY=:1 $1xfoil < $2 > $3
fi
kill -15 $xvfb_pid &> /dev/null
