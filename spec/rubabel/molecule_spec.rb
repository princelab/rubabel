require 'spec_helper'

require 'rubabel/molecule'

describe Rubabel::Molecule do
  before(:each) do
    @mol = Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' )
  end

  attributes = {
    charge: 0,
    formula: "C27H46O",
    dim: 2,
    spin: 1,
  }.each do |methd, exp|
    it "has #{methd}" do
      @mol.send(methd).should == exp
    end
  end

  it '#each iterates through each atom in id order' do
    cnt = 0
    @mol.each do |atom|
      atom.id.should == cnt
      cnt += 1
    end
    @mol.atoms.size.should == 33
    @mol.add_h!
    @mol.atoms.size.should == 74
  end

  it '#hydrogens_added?' do
    @mol.hydrogens_added?.should be_false
    @mol.atoms.size.should == 33
    @mol.add_h!
    @mol.atoms.size.should == 74
    @mol.hydrogens_added?.should be_true
  end

  it 'calculates #ob_sssr (smallest set of smallest rings)' do
    ar = @mol.ob_sssr
    ar.should be_an(Array)
    # in the future should be Rubabel::Ring
    ar.first.should be_a(OpenBabel::OBRing)
  end

  describe 'getting other descriptors' do

  end

  describe '3D' do
    before(:each) do
      # cholesterol
      @mol = Rubabel::Molecule.from_string("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    end

    it 'can be turned into a 3D molecule' do
      # this is only as good as the Builder is good.  For instance, it fails
      # to get all the stereo centers of cholesterol (but it does warn on
      # this, although I don't know how to capture the warnings (can't get
      # with stdout or stderr??))
      @mol.ob.has_3d.should be_false
      @mol.make_3d!
      @mol.ob.has_3d.should be_true 
    end
  end

  describe 'breaking a molecule' do
    before(:each) do
      @mol =  Rubabel::Molecule.from_string("NC(=O)CO")
    end

    it 'num_atoms, atoms and each_atom are sensitive to #add_h!' do
      @mol.num_atoms.should == 5
      @mol.atoms.size.should == 5
      @mol.each_atom.map.to_a.size.should == 5
      @mol.add_h!
      @mol.num_atoms.should == 10
      @mol.atoms.size.should == 10
      @mol.each_atom.map.to_a.size.should == 10
    end

    it 'num_bonds, bonds and each_bond are also sensitive to #add_h!' do
      @mol.num_bonds.should == 4
      @mol.bonds.size.should == 4
      @mol.each_bond.map.to_a.size.should == 4
      @mol.add_h!
      @mol.num_bonds.should == 9
      @mol.bonds.size.should == 9
      @mol.each_bond.map.to_a.size.should == 9
    end

    it 'can be split into multiple molecules' do
      reply = @mol.split(@mol.bonds.first, @mol.bonds.last)
      reply.should be_a(Array)
      reply.size.should == 3
      csmiles = reply.map(&:csmiles)
      csmiles.sort.should == %w(N CC=O O).sort
    end
  end

  describe 'matching patterns (SMARTS)' do
    before(:each) do
      @mol =  Rubabel::Molecule.from_string("NC(=O)O")
    end

    it 'can match smarts patterns' do
      smarts_pattern = 'C=O'
      ar = @mol.matches(smarts_pattern)
      ar.should be_a(Array)
      ar.size.should == 1
      ar.first.map(&:type).should == %w(Cac O.co2)

      # reverse order of smiles gives reverse order
      ar = @mol.matches(smarts_pattern.reverse)
      ar.first.map(&:type).should == %w(Cac O.co2).reverse
    end

  end
end
