# Irreversible package for [BEAST 2](http://www.beast2.org).

Bayesian irreversible models for developmental biology




## Install

Get BEAST 2 from [beast2.org](http://www.beast2.org).

Get the BEASTLabs package: open BEAUti, menu `File => Manage packages`. Select BEASTlabs from the list in the dialog that pops up and click the `install` button.

Build the Irreversible package as follows: This requires java [JDK 8](https://adoptium.net/temurin/releases/?version=8) (higher java versions may work, lower versions don't).
You also need `ant` available from [here](https://ant.apache.org/bindownload.cgi) and `git`.

Get necessary java code:

* `git clone --depth 1 https://github.com/CompEvol/beast2.git`
* `git clone --depth 1 https://github.com/BEAST2-Dev/BEASTLabs.git`
* `git clone https://github.com/rbouckaert/irreversible.git`

Compile with `ant`. 

* `cd BEASTLabs`
* `ant addon`
* `cd ../irreversible`
* `ant addon`

Go to the BEAST package directory  (~/.beast/2.6/ for Linux, ~/Library/Application Support/BEAST/2.6/ on OS X, or c:\Users\yourname\BEAST\2.6\ on Windows). 

* `mkdir irreversible`
* `cd irreversible`
* `cp <path to irreversible>/build/dist/Irreversible.addon.v0.0.1.zip .`
* `unzip Irreversible.addon.v0.0.1.zip`

If you run `<path to beast>/bin/packagemanager -list` you should have `irreversible` as one of the packages listed.

## Setting up analyses with different data

* load alignment in BEAUti, save to XML
* copy/paste all `sequences` elements into example XML
* run BEAST on the edited XML file
