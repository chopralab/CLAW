#!/bin/bash
# Script to launch LipidLynxX to allow the user to upload dataset from different naming conventions and levels

/lipidmaps/python-3.8.4/bin/python3 /lipidmaps/LipidLynxX/cli_lynx.py convert-file $1 --column 0 --output $2 --style BioPAN
