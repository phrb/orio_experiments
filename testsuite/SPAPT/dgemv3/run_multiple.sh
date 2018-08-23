#!/bin/bash
set -e

for i in `seq 1 $1`;
do
	../../../orcc -v -e $2
done
