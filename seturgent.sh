#!/bin/sh

IDS=$(xdotool search --name 'MatlabT')

for id in $IDS; do
    seturgent "${id}"
done
