---
layout: page
title: How to connect to a remote computer
order: 2
session: 1
length: 20
toc: false
---

## Remote Computing

Bioinformatic analysis can require huge amounts of computational resource. It is often impractical to use acquire enough resource to install on or under your desk!

For this reason, we usually connect to a larger server (or cluster of servers) in the University's data-centre, or possibly in the cloud.

Data and applications are stored and processed on remote servers. Users can access these resources and perform tasks remotely using various devices, such as computers, smartphones, or tablets, via an internet connection.

## Command line interface

Much bioinformatic analysis is performed from the command line on a Linux server. But first we need to use a program to connect to the remote server.  This is often refered to as the 'host' and to connect you will need to know its name or IP address.
e.g. for ISCA the hostname is `login.isca.ex.ac.uk`

We will connect using ssh (Secure SHell) On a Mac you can use this from the terminal directly, from a Windows based machine you need to install a program.

## Windows Applicaton style options are

### MobaXterm

[MobaXterm](https://mobaxterm.mobatek.net) has a free portable version.

### Termius

[Termius](https://termius.com/) is available as a free application on Microsoft Store

### PuTTy

[PuTTy](https://www.putty.org/) used to be the 'go-to' option and is still very popular. I can neither recommend nor ignore it!

## Local command line options

### Windows

Computers with Windows operating systems do not automatically have a Unix Shell program
installed.

### MacOS

For a Mac computer running macOS Mojave or earlier releases, the default Unix Shell is Bash.
For a Mac computer running macOS Catalina or later releases, the default Unix Shell is Zsh.
Your default shell is available via the Terminal program within your Utilities folder.

To open Terminal, try one or both of the following:

* In Finder, select the Go menu, then select Utilities.
  Locate Terminal in the Utilities folder and open it.
* Use the Mac 'Spotlight' computer search function.
  Search for: `Terminal` and press Return>.

To check if your machine is set up to use something other than Bash,
type `echo $SHELL` in your terminal window.

If your machine is set up to use something other than Bash,
you can run it by opening a terminal and typing `bash`.

Here are instruction on [how to Use Terminal on a Mac](http://www.macworld.co.uk/feature/mac-software/how-use-terminal-on-mac-3608274/)

## Example

See this page from the RSE workshop with more details on [Connecting to ISCA](https://uniexeterrse.github.io/intro-to-isca/03_connection/index.html)

\* material sourced from [RSE pages](https://uniexeterrse.github.io/intro-unix-shell/setup.html)
