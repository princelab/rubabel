require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/atom'

describe Rubabel::Atom do

  it 'can be created given an element symbol' do
    hydrogen = Rubabel::Atom[:h]
    hydrogen.el.should == :h
    hydrogen.id.should == 0

    carbon = Rubabel::Atom[:c]
    carbon.el.should == :c
    carbon.id.should == 0

    chlorine = Rubabel::Atom[:cl, 3]
    chlorine.el.should == :cl
    chlorine.id.should == 3
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
