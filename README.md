# ClInt


## Installation Procedure

For simple installation, if all preqrequisite programs are already installed, a normal `git clone` suffices. If the need arises to clone /some/ of the needed programs, the command below can be used

`git clone --recurse-submodules --shallow-submodules https://github.com/dmgie/ClInt.git`

Or

`git clone --recurse-submodules --shallow-submodules git@github.com:dmgie/ClInt.git`

This will clone some of the necessary programs into a subfolder, at least the ones available from github. Using the optional `install_subprograms.sh` script, it will try to compile the cloned submodules. For programs without a repo (e.g only a binary available) it will try to install these using `wget`. 

Other (such as ones listed directly below this line) might require some manual intervention).
- Example Program one
- Example Program two

## How To Use It

The workflow file can be called using the `workflow.bash` script which allows selection programs to run for the pipeline. 

<Example Image>

## Programs Required:
- `minimap2`
- `fastp`
- `samtools`
- `bcftools`
- `MUSIAL`
