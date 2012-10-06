require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/atom'

describe Rubabel::Atom do

  it 'can be created given an element symbol' do
    hydrogen = Rubabel::Atom[:h]
    hydrogen.el.should == :h

    carbon = Rubabel::Atom[:c]
    carbon.el.should == :c

    chlorine = Rubabel::Atom[:cl, 3]
    chlorine.el.should == :cl
    chlorine.id.should == 3
  end

  specify 'equality stuff!!'

  it 'properly removes hydrogens (explicit H)' do
    mol = Rubabel["CC"]
    mol.add_h!
    mol.atoms[0].remove_an_h!
    mol.atoms[0].charge += 1
    mol.formula.should == "C2H5"
    mol.csmiles.should == 'C[CH2+]'
    mol.exact_mass.round(5).should == 29.03913
    mol.charge.should == 1

    # need to get adding working!!!
    #p mol.atoms[0].spin
    #mol.atoms[0].add_an_h!
    #p mol.atoms[0].spin
    #mol.formula.should == 'C2H6'
    #mol.atoms[0].charge -= 1
    #mol.exact_mass.should == 323
    #mol.charge.should == 0
  end

  it 'can find atom identities with simple questions' do
    mol = Rubabel["NCC(O)CC(=O)"]
    (c_exp, o_exp) = mol.matches("C=O").first
    mol.find(&:carbonyl_carbon?).id.should == c_exp.id
    mol.find(&:carbonyl_oxygen?).id.should == o_exp.id
  end

  it 'properly removes and adds hydrogens (implicit H)' do
    mol = Rubabel["CC"]
    mol.atoms[0].remove_an_h!
    mol.formula.should == "C2H5"
    mol.csmiles.should == 'C[CH2]'
    mol.atoms[0].charge += 1
    mol.csmiles.should == 'C[CH2+]'

    # still need to properly implement
    #mol.atoms[0].add_an_h!
    #mol.formula.should == 'C2H6'
    #p mol
  end

  describe 'working with a complex molecule' do

    before do
      @mol = Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' )
      @atom = @mol.atoms.first
      @mol_h =  Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' )
      @mol_h.add_h!
      @atom_with_h = @mol_h.atoms.first
    end

    attributes = {
      charge: 0,
      id: 0,
      spin: 0,
    }.each do |methd, exp|
      it "has #{methd}" do
        @atom.send(methd).should == exp
      end
    end


    it '#mol retrieves the parent molecule' do
      @atom.mol.should == @mol

      # no parent molecule 
      h = Rubabel::Atom[:h]
      h.mol.should be_nil
    end

    it 'can get the bonds' do
      @atom.each_bond do |bond|
        bond.should be_a(Rubabel::Bond)
      end
      @atom.bonds.size.should == 4
    end

    it 'can add a bond' do
    end

    it 'can get the neighboring atoms' do
      @atom.id.should == 0
      @atom.atomic_num.should == 6
      @atom.type.should == 'C3'
      @atom.each_atom do |nbr_atom|
        nbr_atom.should be_a(Rubabel::Atom)
        nbr_atom.ob.equal?(@atom.ob).should be_false
      end
      @atom.atoms.size.should == 4
    end

    it '#get_bond can retrieve a particular bond based on the atom' do
      other = @mol.atoms[3]
      bond = @atom.get_bond(other)
      bond.atoms.map(&:id).should == [0,3]
      other = @mol.atoms[15]  # these are not connected
      @atom.get_bond(other).should be_nil
    end

    it '#coords gets the coordinates' do
      @atom.coords
    end
  end
end
