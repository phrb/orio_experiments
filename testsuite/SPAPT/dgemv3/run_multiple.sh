#!/bin/bash
for i in `seq 1 $1`;
do
	../../../orcc -v -e atax2.src6.c
done    
