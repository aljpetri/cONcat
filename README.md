# cONcat: Computational reconstruction of concatenated fragments from long Oxford Nanopore reads

## Installation Guide <a name="installationguide"></a>
At the moment building from source is the only option to install the tool. This requires users to install the Rust programming language onto their system.

## Installing Rust <a name="installingrust"></a>
You can install rust via<br />

`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh` (for macOS and Linux or other Unix-based OS). For Windows please follow the instructions on the following site: https://forge.rust-lang.org/infra/other-installation-methods.html .<br />

## Installation <a name="installation"></a>
After cloning the repository via `git@github.com:aljpetri/cONcat.git` use the following two commands to compile the code: <br />
`cd cONcat/cONcat_Code` <br />
`cargo build --release` ( Compile the current package, the executable is then located in target/release) <br />
Please note that cONcat depends on [Edlib](https://github.com/Martinsos/edlib) and, therefore, has dependencies to cmake and g++ . 


## Running
`cONcat  --expected Path/to/Expected_fragments.csv --fastq path/to/1.fastq --outfile path/to/outfile --verbose`


## Credits <a name="credits"></a>

Please cite this study when using cONcat: 

cONcat: Computational reconstruction of concatenated fragments from long Oxford Nanopore reads, 
Alexander J. Petri, Mai Thi-Huyen Nguyen, Anjali Rajwar, Erik Benson, Kristoffer Sahlin [LINK](https://www.biorxiv.org/content/10.1101/2025.03.05.641699v2.abstract).
