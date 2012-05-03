require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/atom'

describe Rubabel::Atom do
  subject { Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' ).atoms.first }

  attributes = {
    charge: 0,
    id: 0,
    spin: 0,
  }.each do |methd, exp|
    it "has #{methd}" do
      subject.send(methd).should == exp
    end
  end

  it 'can get the bonds' do
    atom = subject
    atom.each_bond do |bond|
      bond.should be_a(Rubabel::Bond)
    end
    # is this correct???
    atom.bonds.size.should == 3
    warn 'warn: not sure if correct'
  end

  it 'can get the neighboring atoms' do
    atom = subject
    atom.each_atom do |nbr_atom|
      nbr_atom.should be_a(Rubabel::Atom)
      nbr_atom.obatom.equal?(atom.obatom).should be_false
    end
    # is this correct???
    atom.atoms.size.should == 4
    warn 'warn: not sure if correct'
  end
end
