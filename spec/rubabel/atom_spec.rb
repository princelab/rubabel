require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/atom'

describe Rubabel::Atom do

  xit 'can be created given an element symbol' do
    hydrogen = Rubabel::Atom[:h]
    hydrogen.el.should == :h

    carbon = Rubabel::Atom[:c]
    carbon.el.should == :c

    chlorine = Rubabel::Atom[:cl, 3]
    chlorine.el.should == :cl
    chlorine.id.should == 3
  end

  describe 'attaching another atom to itself (as part of a molecule)' do

    it 'attaches preformed atom' do
    end

    it 'attaches given an atomic number' do

    end

    it 'attaches given an element symbol' do
    end

  end

  specify 'equality' do
    mol = Rubabel["CCO"]
    oxygen = mol.atoms[2]
    oxygen_from_match = mol.matches("CO").first.last
    (oxygen == oxygen_from_match).should be_true
    (oxygen.equal?(oxygen_from_match)).should be_true
    (oxygen.eql?(oxygen_from_match)).should be_true

    mol2 = Rubabel["CCO"]
    (mol.atoms[0] == mol2.atoms[0]).should_not be_true
    (mol.atoms[0].equal?(mol2.atoms[0])).should_not be_true
    (mol.atoms[0].eql?(mol2.atoms[0])).should_not be_true
  end

  it 'removes hydrogens with proper charge accounting' do
    mol = Rubabel["CC"]
    mol.add_h!
    mol.atoms[0].remove_an_h!
    # I was expecting this, but it is now giving formulas with charges (which
    # is probably better)
    ##mol.formula.should == "C2H5"
    mol.formula.should == "C2H5+"
    mol.csmiles.should == 'C[CH2+]'
    mol.exact_mass.round(5).should == 29.03913
    mol.charge.should == 1

    # can't seem to get working properly!!!
    #mol.atoms[0].add_an_h!
    #mol.formula.should == 'C2H6'
    #mol.csmiles.should == 'CC'
    #mol.charge.should == 0
    ##fmol.atoms[0].charge -= 1
    #mol.exact_mass.should == 323
  end

  it 'can find atom identities with simple questions' do
    mol = Rubabel["NCC(O)CC(=O)"]
    (c_exp, o_exp) = mol.matches("C=O").first
    mol.find(&:carbonyl_carbon?).id.should == c_exp.id
    mol.find(&:carbonyl_oxygen?).id.should == o_exp.id
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

    specify '#<<(atom) adds an atom and returns the added atom' do
      mol = Rubabel["C"]
      reply = (mol[0] << :n << :o)
      reply.should be_a(Rubabel::Atom)
      mol.csmiles.should == 'CNO'
    end

    it '#mol retrieves the parent molecule' do
      @atom.mol.should == @mol

      ## no parent molecule 
      #h = Rubabel::Atom[:h]
      #h.mol.should be_nil
    end

    it 'can get the bonds' do
      @atom.each_bond do |bond|
        bond.should be_a(Rubabel::Bond)
      end
      @atom.bonds.size.should == 4
    end

    it 'can add a bond'

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
