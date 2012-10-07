# rubabel

Ruby interface to the openbabel ruby bindings (or the openbabel gem).  The
interface attempts to be a ruby-ish analogue of
[pybel](http://openbabel.org/docs/current/UseTheLibrary/Python_PybelAPI.html).

## Examples

The [Chemistry Toolkit Rosetta Wiki](http://ctr.wikia.com/wiki/Chemistry_Toolkit_Rosetta_Wiki) has a lot of examples you can check out.

### Creating a molecule

#### From a string

    require 'rubabel'

    # by default, reads in smiles strings
    serine = Rubabel["C(C(C(=O)O)N)O"]
    # more formally:
    serine = Rubabel::Molecule.from_string("C(C(C(=O)O)N)O")

    # also any other format openbabel supports, for example inchi
    serine = Rubabel["InChI=1S/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)", :inchi]

    # from the internet:
    mol = Rubabel[some_molecule, Rubabel.format_from_mime(some_mime_type)]

Find out all the formats Rubabel supports (hash is format key pointing to the description):

    hash = Rubabel.in_formats
    hash = Rubabel.out_formats

#### From a file

Reading multiple entries from a file:
    
    Rubabel.foreach("file.sdf") do |mol|
      puts mol.exact_mass
    end

Foreach returns an enumerator that can be chained:

    # return an array of every unique atom type in the file
    uniq_atom_types = Rubabel.foreach("file.mol").flat_map {|mol| mol.map(&:type) }.uniq

Read a single molecule from a file (reads only the first molecule)

    mol = Rubabel::Molecule.from_file("file.sdf")
    # handles gzipped files seamlessly:
    mol = Rubabel::Molecule.from_file("file.sdf.gz") 
    mol = Rubabel.molecule_from_file("file.sdf") # alternative

    # explicit format for uninformative/wrong extensions:
    mol = Rubabel::Molecule.from_file("file", :sdf)

### Writing/Drawing

#### create string output

    mol = Rubabel["OCC"] # ethanol

    mol.to_s       # canonical smiles -> "CCO"
    mol.csmiles    # same thing

    mol.to_s(:smi) # smiles -> "OCC"
    mol.smiles     # same thing

For inclusion in a file with standard smiles formatting (SMILES\tID\n):

    can_smiles_string = mol.write_string # -> "CCO\t\n"
    mol.title = "ethanol"
    can_smiles_string = mol.write(:can)  # -> "CCO\tethanol\n"

Other formats in the same manner:

    pdb_string = mol.write(:pdb)

Write to a file directly (single molecule only; depends on file extension for type):

    # write to a smiles file
    mol.write("file.smi")
    mol.write_file("file.smi")

Write multiple molecules to a file:

    File.open("somefile.pdb", 'w') do |out|
      molecules.each {|mol| out.print mol.write(:pdb) }
    end

#### Drawing

If you write to svg or png (png uses mini_magick to convert from svg) then the
molecule is automatically drawn in 2D:

    mol = Rubabel["NCC(O)C(=O)O"]
    mol.write("file.svg")

    # must have imagemagick ('convert' command) and mini_magick gem installed
    mol.write("file.png") 

### Searching and Splitting

*each_match*, *matches*, *matches?*, *smarts_indices* all take the same input (SMARTS
string or object and optional boolean specifying uniqueness of results):

    mol = Rubabel["NCC(O)C(=O)O"]
    mol.each_match("CO") do |match|
      # match is just an array of atoms that matched
      match.first.el # => :c
      match.last.el # => :o
    end
    
    # matches returns all the matches in an array
    all_matches = mol.matches("CO")
    # all the match routines take a boolean to alter uniqueness
    all_matches = mol.matches("CO", false) # some matches may not be uniq

Have some bonds to break?, split makes new molecules split from that bond(s)
    
    bonds = mol.matches("CO").map {|c, o| c.get_bond(o) }
    mol.split(*bonds)  # splits between every carbon single bonded to oxygen

### Add & Delete atoms/bonds

    mol = Rubabel["OCC"]
    carbon = mol.add_atom!(6)
    # add the new carbon to the last carbon
    last_carbon = mol[2]
    # add it directly to the carbon:
    last_carbon.add_atom!(carbon)
    # alternatively, add it at the molecule level:
    mol.add_bond!(last_carbon, carbon, 2)


    


## Installing

First, many thanks to Andreas Maunz for packaging openbabel as a gem which makes this install quite painless.

### Quick Install

On a POSIX system, make sure you have openbabel (including header files), cmake, curl, tar, sed and make {see openbabel instructions}[https://github.com/amaunz/openbabel-gem].  On ubuntu/debian:

    sudo apt-get install openbabel libopenbabel-dev cmake make curl

Then install the gem (which should install the openbabel gem, too):

    gem install rubabel

### Building from Source

1. download openbabel 
2. swap out Init_OpenBabel for Init_openbabel in scripts/ruby/openbabel-ruby.cpp (see [here](http://forums.openbabel.org/Ruby-Open-Babel-in-2-1-1-td957640.html)).  Some versions have this fixed already, apparently.
3. make sure you have the right [dependencies to compile](http://openbabel.org/docs/2.3.1/Installation/install.html#compiling-open-babel)

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
