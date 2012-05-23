require 'spec_helper'

require 'rubabel'

describe Rubabel::Molecule::Fragmentable do
  describe 'the :co rule' do

    describe 'water loss' do

      it ':h2oloss' do
        mol = Rubabel["NCCC(O)CC"]
        fragments = mol.fragment( rules: [:h2oloss] )
        fragments.map(&:csmiles).sort.should == ["CC=CCC[NH3+]", "CCC=CC[NH3+]", "O", "O"]
      end

      it ':h2oloss [does not allow bad chemistry]' do
        # lone pair and double bond resonance ?
        mol = Rubabel["NCC(O)CC"]
        fragments = mol.fragment( rules: [:h2oloss] )
        fragments.map(&:csmiles).sort.should == ["CC=CC[NH3+]", "O"]

        mol = Rubabel["NC(O)CC"] 
        fragments = mol.fragment( rules: [:h2oloss] )
        fragments.map(&:csmiles).sort.should == []
      end
    end

    describe 'backbone cleavage' do

      it 'cleaves beside alcohols yielding ' do
        mol = Rubabel["NCCC(O)CC"]
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
        mol = Rubabel["NCCC(=O)OC"]
        pieces = mol.fragment( rules: [:co] )
        pieces.should be_empty
      end
    end
  end
end
