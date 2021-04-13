#!/bin/bash
MY_DIRNAME=$(dirname $0)
cd $MY_DIRNAME
find . -name "main" -exec rm -r {} \;
find . -name "main.o" -exec rm -r {} \;
make && ./main