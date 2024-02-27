#!/bin/bash

#Script to push all added scripts to this directory

git add .
git commit -m "push-scripts"
git branch -M main
git push -u origin main
