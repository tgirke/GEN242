#!/bin/bash
## Convenience script to add, commit and push updates to GitHub 
## Author: Thomas Girke
## Last update: Mar 21, 2024
## Note: for creating a new instance of Hugo site prior to quarter, do the following:
## (1) create a new (empty) repos on GitHub named GEN242-<YEAR>; (2) Clone repos to local; (3) Copy 
## content of latest GEN242 version into new instance with: "cp -rT ../GEN242 ."; (4) Update
## .git/config with correct path here GEN242-<YEAR>.git. (5) Run git add -A :/; git commit -am "..."; git push
## (6) Enable online GitHub Pages under Settings/Pages where you use the settings from original instance (here GEN242) 

## (1) Makes sure you are in correct branch (currently there is only one) and remote changes are pulled
git checkout main  
git pull

## (2) Build site and copy rendered pages to docs 
# rm -rf public/* # clean out public -> not required with publishDir = "docs" in config.toml
# rm -rf docs/* # clean out docs (assumes dir exists) -> not required with publishDir = "docs" in config.toml
# Rscript -e "blogdown::build_site()" # Build site with blogdown via command-line call
Rscript -e "blogdown::build_site(build_rmd = 'md5sum')" # runs build only for rmd pages that have changed
# cp -r public/* docs/ # Copy build/rendered pages to docs -> not required with publishDir = "docs" in config.toml

## (3) Commit and push changes
git add -A :/
git commit -am "some edits"
git push -u origin main
echo "Committed/pushed changes to GitHub"
