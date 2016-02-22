#!/bin/bash
# 
# Convert all pdf figures to png for convenience. Requires convert 
# from ImageMagick
# 

FIGDIR="$1"

if [ -z "$FIGDIR" ]; then 
  echo "Provide the figures directory"
  exit 1
fi

# Go ahead
for file in $(find "$FIGDIR" -name "*.pdf"); do 
  
  echo $file
  convert -density 180 "$file" -flatten "${file%%.pdf}.png"
  
done
