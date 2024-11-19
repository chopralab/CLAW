#!/bin/bash

# Script Name: unstage_png_batches.sh
# Description: Unstages .png files in batches of 10 to prevent system overload.
# Usage: ./unstage_png_batches.sh

echo "---------------------------------------------"
echo "Starting the unstage_png_batches.sh script..."
echo "---------------------------------------------"

# Initialize an empty array to hold filenames
files=()

# List all staged .png files with null separators to handle special characters
echo "Listing all staged .png files..."
git ls-files --cached -- '*.png' -z | while IFS= read -r -d '' file; do
    files+=("$file")
    
    # When 10 files are collected, unstage them
    if [ ${#files[@]} -eq 10 ]; then
        echo "---------------------------------------------"
        echo "Unstaging the following 10 files:"
        for f in "${files[@]}"; do
            echo "  $f"
        done
        echo "---------------------------------------------"
        
        git reset HEAD -- "${files[@]}"
        
        if [ $? -eq 0 ]; then
            echo "Successfully unstaged 10 files."
        else
            echo "Error: Failed to unstage some files."
        fi
        
        # Reset the array
        files=()
    fi
done

# Unstage any remaining files (less than 10)
if [ ${#files[@]} -gt 0 ]; then
    echo "---------------------------------------------"
    echo "Unstaging the remaining ${#files[@]} file(s):"
    for f in "${files[@]}"; do
        echo "  $f"
    done
    echo "---------------------------------------------"
    
    git reset HEAD -- "${files[@]}"
    
    if [ $? -eq 0 ]; then
        echo "Successfully unstaged ${#files[@]} file(s)."
    else
        echo "Error: Failed to unstage some files."
    fi
fi

echo "---------------------------------------------"
echo "Completed unstaging .png files."
echo "---------------------------------------------"

# Optional: Count how many .png files are now unstaged
unstaged_count=$(git ls-files --others --exclude-standard -- '*.png' | wc -l)
echo "Total .png files now untracked: $unstaged_count"
