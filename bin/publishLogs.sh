#!/bin/bash
ID=$1

mkdir -p logs
cp .command.log  logs/${ID}.log  
cp .command.err  logs/${ID}.err
cp .command.out  logs/${ID}.out
cp .command.run  logs/${ID}.run
