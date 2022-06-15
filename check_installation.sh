#!/bin/bash

# Check if fastp is installed
conda list -n wgregseq fastp

# Check if bbmap is in the correct location
FILE=./bbmap/bbmap.sh
if test -f "$FILE"; then
    echo "$FILE exists."
else
    read -p "$FILE does not exist! Want to download bbmap?.(Yy/Nn)" yn
    case $yn in
        [Yy]* ) {
            wget 'https://sourceforge.net/projects/bbmap/files/latest/download/BBMap_38.96.tar.gz';tar -xf 'BBMap_38.96.tar.gz' ;rm 'BBMap_38.96.tar.gz';
        } || {
            echo "installing bbmap failed, please download it manually and store it in this repository."
            };;
        [Nn]* ) echo "BBmap was not installed. Please install manually.";;
        * ) echo "Please answer yes or no.";;
    esac 
fi

{ # try
    dirname $(greadlink -f $0) &&
    echo "greadlink works"

} || { 
    read -p "greadlink does not work! Do you want to install it.(Yy/Nn)" yn
    case $yn in
        [Yy]* ) {
            brew install 'coreutils' && echo "greadlink should work now";
        } || {
            read -p "brew is not installed! Do you want to install it? (Yy/Nn)" yn
            case $yn in
                [Yy]* ) /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)";brew install 'coreutils' && echo "greadlink should work now";;
                [Nn]* ) echo "Not installing brew. Please make sure that greadlink works manually."; break;;
                * ) echo "Please answer yes or no";;
            esac
        };;
        [Nn]* ) echo "Please make sure that greadlink works. Can be done by using `brew coreutils`";;
        * ) echo "Please answer yes or no.";;
    esac
}