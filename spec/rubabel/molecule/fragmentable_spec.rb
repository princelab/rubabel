require 'spec_helper'

require 'rubabel'

describe Rubabel::Molecule::Fragmentable do
  describe 'the :co rule' do

    describe 'water loss' do

      it ':h2oloss' do
        mol = Rubabel["NCCC(O)CC"]
        fragments = mol.fragment( rules: [:h2oloss] )
        fragments.flatten(1).map(&:csmiles).sort.should == ["CC=CCC[NH3+]", "CCC=CC[NH3+]", "O", "O"]
      end

      it ':h2oloss [does not allow bad chemistry]' do
        # lone pair and double bond resonance ?
        mol = Rubabel["NCC(O)CC"]
        fragments = mol.fragment( rules: [:h2oloss] )
        fragments.flatten(1).map(&:csmiles).sort.should == ["CC=CC[NH3+]", "O"]

        mol = Rubabel["NC(O)CC"] 
        fragments = mol.fragment( rules: [:h2oloss] )
        fragments.flatten(1).map(&:csmiles).sort.should == []
      end
    end

    describe 'backbone cleavage' do

      it 'cleaves beside alcohols yielding aldehydes' do
        mol = Rubabel["NCCC(O)CC"]
        mol.correct_for_ph!
        total_mass = mol.add_h!.mass

        pieces = mol.fragment(rules: [:co])
        pieces.size.should == 2
        pieces.map(&:size).should == [2,2]
        pieces.flatten(1).map(&:csmiles).should == ["CC[NH3+]", "CCC=O", "C(C=O)C[NH3+]", "CC"]
        pieces.each do |pair|
          pair.map(&:mass).reduce(:+).should == total_mass
        end
      end

      xit 'does not cleave esters' do
        mol = Rubabel["NCCC(=O)OC"]
        pieces = mol.fragment( rules: [:co] )
        pieces.should be_empty
      end
    end
  end
end
