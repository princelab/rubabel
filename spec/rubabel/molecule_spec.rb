require 'spec_helper'

require 'rubabel/molecule'

describe Rubabel::Molecule do
  subject { Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' ) }

  attributes = {
    charge: 0,
    formula: "C27H46O",
    dim: 2,
    spin: 1,
  }.each do |methd, exp|
    it "has #{methd}" do
      subject.send(methd).should == exp
    end
  end

  it '#each iterates through each atom in id order' do
    cnt = 0
    subject.each do |atom|
      atom.id.should == cnt
      cnt += 1
    end
    subject.atoms.size.should == 33
    subject.add_h!
    subject.atoms.size.should == 74
  end

  it 'calculates #ob_sssr (smallest set of smallest rings)' do
    ar = subject.ob_sssr
    ar.should be_an(Array)
    # in the future should be Rubabel::Ring
    ar.first.should be_a(OpenBabel::OBRing)
  end

  describe 'breaking a molecule' do
    subject { Rubabel::Molecule.from_string("NC(=O)O" ) }

    it 'num_atoms, atoms and each_atom are sensitive to #add_h!' do
      subject.num_atoms.should == 4
      subject.atoms.size.should == 4
      subject.each_atom.map.to_a.size.should == 4
      subject.add_h!
      subject.num_atoms.should == 7
      subject.atoms.size.should == 7
      subject.each_atom.map.to_a.size.should == 7
    end

    it 'num_bonds, bonds and each_bond are also sensitive to #add_h!' do
      subject.num_bonds.should == 3
      subject.bonds.size.should == 3
      subject.each_bond.map.to_a.size.should == 3
      subject.add_h!
      subject.num_bonds.should == 6
      subject.bonds.size.should == 6
      subject.each_bond.map.to_a.size.should == 6
    end

    it 'can be split into smaller molecules' do
      reply = subject.split(subject.bonds.first)
      reply.should be_a(Array)
      reply.size.should == 2
      reply.any? {|mol| mol.csmiles == 'N'}.should be_true
      reply.any? {|mol| mol.csmiles == 'OC=O'}.should be_true
    end

    it 'can match smarts patterns' do
      puts "HERE"
      p subject.matches('C=O')
      p subject.matches('O=C')
    end

  end
end
