# rubabel

Ruby interface to the openbabel ruby bindings (or the openbabel gem).  The
interface attempts to be a ruby-ish analogue of
[pybel](http://openbabel.org/docs/current/UseTheLibrary/Python_PybelAPI.html).

## Examples

The [Chemistry Toolkit Rosetta Wiki](http://ctr.wikia.com/wiki/Chemistry_Toolkit_Rosetta_Wiki) has a lot of examples you can check out.

### Creating a molecule

    require 'rubabel'

    # by default, reads in smiles strings
    serine = Rubabel["C(C(C(=O)O)N)O"]

    # also any other format openbabel supports (a ton), for example inchi
    serine = Rubabel["InChI=1S/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)", :inchi]

    
    Rubabel.foreach("file.sdf") do |mol|
      print mol.write  # -> canonical smiles: "smiles_string\tname\n"
      puts mol  # .to_s -> canonical smiles string with no name or newline
      puts mol.exact_mass
    end

    uniq_atom_types = Rubabel.foreach("file.mol").map(&:type).uniq

## Installing

First, many thanks to Andreas Maunz for packaging openbabel as a gem which makes this install quite painless.

### Quick Install

On a POSIX system, make sure you have openbabel (including header files), cmake, curl, tar, sed and make {see openbabel instructions}[https://github.com/amaunz/openbabel-gem].  On ubuntu/debian:

    sudo apt-get install openbabel libopenbabel-dev cmake make curl

Then install the gem (which should install the openbabel gem, too):

    gem install rubabel

### Building from Source

1. download openbabel 
2. swap out Init_OpenBabel for Init_openbabel in scripts/ruby/openbabel-ruby.cpp (see here[http://forums.openbabel.org/Ruby-Open-Babel-in-2-1-1-td957640.html]).  Some versions have this fixed already, apparently.
3. make sure you have the right {dependencies to compile}(http://openbabel.org/docs/2.3.1/Installation/install.html#compiling-open-babel)

Here's a complete example of compiling for a single user on Ubuntu 11.10 and probably will be generally forward compatible for some time.  This will compile bindings on whichever ruby comes up with '/usr/bin/env ruby':

    # install the dependencies:
    sudo apt-get install libeigen2-dev cmake libwxgtk2.8-dev libxml2-dev libcairo2-dev
    # unpack it:
    tar -xzvf openbabel-2.3.1.tar.gz
    # swap out buggy lines in ruby bindings:
    sed -i 's/Init_OpenBabel/Init_openbabel/g' openbabel-2.3.1/scripts/ruby/openbabel-ruby.cpp
    # make a separate build directory for building in:
    mkdir build-rvmruby1.9.3
    cd build-rvmruby1.9.3
    mkdir ~/tools
    cmake ../openbabel-2.3.1 -DRUBY_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/tools/openbabel-rvmruby1.9.3
    make && make install

[[Still need directions to install the gem on top of a build from source]]
 
## Copyright

MIT License.  See LICENSE for further details.
