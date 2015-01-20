#!/bin/bash

rm rd2_view
make -f makerd
rm dynQuad
make -f makeTestQuad
./dynQuad | ./rd2_view -
