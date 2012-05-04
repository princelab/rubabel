require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/atom'

describe Rubabel::Atom do
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

    it 'can get the bonds' do
      @atom.each_bond do |bond|
        bond.should be_a(Rubabel::Bond)
      end
      @atom.bonds.size.should == 3
    end

    it 'can get the neighboring atoms' do
      @atom.id.should == 0
      @atom.atomic_num.should == 6
      @atom.type.should == 'C3'
      @atom.each_atom do |nbr_atom|
        nbr_atom.should be_a(Rubabel::Atom)
        nbr_atom.obatom.equal?(@atom.obatom).should be_false
      end
      @atom.atoms.size.should == 4
    end

    it '#coords gets the coordinates' do
      @atom.coords
    end
  end
end
