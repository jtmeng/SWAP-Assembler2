#!/bin/bash
g++ -std=c++0x -O3 simulator.cpp -o simu

for ((i=0; i<16; ++i))  
do  
	./simu -mutate 0.05 -err 0.01 -orgRef ../human_ref_withoutM.fa -mutateRef org.fa -errReads $i.fa
done  
