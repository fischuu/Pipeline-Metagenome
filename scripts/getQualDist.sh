#!/bin/sh

zcat $1 | head -n 100000 | awk 'NR == 0 || NR % 4 == 0' | sed 's/\(.\)/\1\n/g' | sort | uniq -ic