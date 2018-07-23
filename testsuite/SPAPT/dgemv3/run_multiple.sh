#!/bin/bash
for i in `seq 1 $1`;
do
	../../../orcc -v -e $2
done    
