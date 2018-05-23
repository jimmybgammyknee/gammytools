#!/bin/bash -l

# check all gzip files in a directory

for i in *.gz
 do 
	echo $i
	gzip -t $i && echo ok || echo bad
done
