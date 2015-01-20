#!/bin/bash

rm rd2_view testBezier
make -f makerd
make -f maketest
./rd2_view s55.rd
./rd2_view s56.rd
./testBezier | ./rd2_view -