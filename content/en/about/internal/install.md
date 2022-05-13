---
title: "Deployment and Maintenance of this Site"
linkTitle: "Installation"
weight: 6
type: docs
---

<br/>
<br/>

> This page provides instructions how to create new deployment instances of this teaching site, and how to configure and customize it. 
> It uses the Docsy theme of the Hugo framework for building content driven websites. 

## Quick start

* Visit this template [page](https://www.dcassol.com/docsy_Rmarkdown/)
* Click on the **Use this Template** button.
* Choose a Repository Name
* Click on the **Create repository from template** button.

![](https://raw.githubusercontent.com/dcassol/images/main/usetemplte.gif)

### Usage locally

* Go to your new repository on GitHub
* Click on the **Code** button.
* Copy the URL `git@github.com:<username>/<repository_name>.git`
* Open terminal and run the following commands.
* Note: every time the repos is cloned from GitHub it is important to include the `--recurse-submodules` argument. This assures that required submodules will be cloned as well. 

```sh
git clone --recurse-submodules --depth 1 git@github.com:<username>/<repository_name>.git
cd <repository_name>
```

* Run the website locally

```sh
hugo server
```

* Run the website locally with `blogdown`

```r
blogdown::serve_site()
```

## Prerequisites and Installation

### Install nodejs

Download `nodejs` binary for 64 bit linux from [here](https://nodejs.org/en/download/). 
Next, install it according the following [instructions](https://bit.ly/3jVJzmU). 
Note, the version in all commands needs to match the downloaded one (here v16.14.2).

1. Unzip the binary archive to any directory, where you wish to install nodejs. The
following uses `/usr/local/lib/nodejs`.

```sh
VERSION=v16.14.2                                                                                                                                                                    
DISTRO=linux-x64                                                                                                                                                                    
sudo mkdir -p /usr/local/lib/nodejs                                                                                                                                                 
sudo tar -xJvf node-$VERSION-$DISTRO.tar.xz -C /usr/local/lib/nodejs
```

2. Set the required environment variables by adding the following lines to your `~/.profile` file.

```sh
# Nodejs
VERSION=v16.14.2                                                                                                                                                                   
DISTRO=linux-x64                                                                                                                                                                    
export PATH=/usr/local/lib/nodejs/node-$VERSION-$DISTRO/bin:$PATH
```

3. Refresh `~/.profile` and test versions

```sh
. ~/.profile
node -v
npm -v
```

4. Make executables available to root when using `sudo npm install`.

```sh
sudo ln -s -f  /usr/local/lib/nodejs/node-v16.14.2-linux-x64/bin/node /usr/local/bin/node                                                                                           
sudo ln -s -f  /usr/local/lib/nodejs/node-v16.14.2-linux-x64/bin/npm /usr/local/bin/npm                                                                                             
sudo ln -s -f  /usr/local/lib/nodejs/node-v16.14.2-linux-x64/bin/npx /usr/local/bin/npx
```

### Install blogdown and Hugo

#### blogdown

```r
installed.packages("rstudio/blogdown")
# If anything wrong try develop version
remotes::install_github("rstudio/blogdown")
```
#### Hugo

You need a recent extended version (we recommend version 0.79.0 or later) of Hugo 
to do local builds and previews of sites that use Docsy.

It is recommended to install `Hugo` from R for working with blogdown

```r
blogdown::install_hugo(extended = TRUE)
```
or from command-line from [here](https://github.com/gohugoio/hugo/releases/)

```bash
wget https://github.com/gohugoio/hugo/releases/download/v0.94.2/hugo_extended_0.94.2_Linux-64bit.deb                                                                                
sudo dpkg -i hugo_extended_0.94.2_Linux-64bit.deb                                                                                                                                   
hugo version
```

For `Windows` and `macOS` please see instructions [here](https://www.docsy.dev/docs/getting-started/). 

#### Install PostCSS


To build or update your site's CSS resources, you also need `PostCSS` to create the final assets. If you need to install it, you must have a recent version of `NodeJS` installed on your machine so you can use `npm`, the Node package manager. By default `npm` installs tools under the directory where you run `npm` install:

```bash
sudo npm install -D autoprefixer
sudo npm install -D postcss-cli
# Starting in version 8 of postcss-cli, you must also separately install postcss:

sudo npm install -D postcss

# go to your website directory
cd <repository_name>
npm audit fix
```

### Run the website locally 

#### with blogdown

* Open R in console or Rstudio

This repo contains an `.Rprofile` file that will automatically serve the site
for you R starting directory is this newly cloned repository. Otherwise, 
change working directory to the repository and run:

```r
blogdown::serve_site()
```

You should see a website is opened in your local browser or Rstudio viewer.

#### Run the website locally from terminal

```bash
cd YOUR_NEW_REPO_PATH
hugo server
```

### Customize

#### Color settings

The background color of the top menu bar can be changed under `assets/scss/_variables_project.scss`. One custom setting in this
file is that the original `$primary: #30638E;` has been changed to `$primary: #28498f;`.

