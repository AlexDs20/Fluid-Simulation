#!/bin/bash

make clean
make
for i in {1..1}
do
  make run
done
