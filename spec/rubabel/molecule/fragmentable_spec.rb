require 'spec_helper'

require 'rubabel'

describe Rubabel::Molecule::Fragmentable do
  describe 'the :co rule' do
    subject { Rubabel::Molecule.from_string("NC(O)CC(=O)O") }
    it 'generates 4 fragments' do
      pieces = subject.fragment(rules: [:co])
      pieces.size.should == 4
      p pieces

      pieces.map(&:csmiles).sort.should == ["O=CCC(=O)[O-]", "O=C[NH2+]", "[NH4+]", "[O-]C(=O)C"]
      pieces.map(&:mass).sort.zip( [18.03382553308, 46.02874015264, 1,1 ]) do |act, exp|
        act.should be_within(0.000001).of(exp)
      end
      pieces.select {|piece| piece.charge > 0 }.size.should == 2
      pieces.select {|piece| piece.charge < 0 }.size.should == 2
    end

    it 'does not cleave esters' do
      mol = Rubabel::Molecule.from_string("NCCC(=O)OC")
      pieces = mol.fragment( rules: [:co] )
      pieces.should be_empty
    end
  end
end
