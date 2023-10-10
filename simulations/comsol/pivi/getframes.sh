#!/bin/bash

for i in *.mp4
do
    name=$(basename "$i" ".mp4")
    mkdir "$name"
    cd "$name"
    ffmpeg -i "../$name.mp4" -r 15 frame%d.png
    cd ..
done
