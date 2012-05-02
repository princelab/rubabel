require 'spec_helper'

require 'rubabel/atom'

describe Rubabel::Atom do
  subject { Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' ).atoms.first }

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

  it 'num_atoms, atoms and each_atom are sensitive to #add_h!' do
    subject.num_atoms.should == 33
    subject.add_h!
    subject.num_atoms.should == 74
  end
  it 'calculates #sssr (smallest set of smallest rings)' do
    ar = subject.sssr
    ar.should be_an(Array)
    # in the future should be Rubabel::Ring
    ar.first.should be_a(OpenBabel::OBRing)
  end
end
