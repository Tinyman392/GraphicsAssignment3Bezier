#!/bin/bash

rm rd2_view
make -f makerd
./rd2_view sqsphere.rd > sphere
#./3view sphere | ./rd_view -
./rd2_view sqtorus.rd > torus
#./3view torus | ./rd_view -
