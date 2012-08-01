require 'spec_helper'

require 'rubabel'

describe Rubabel::Molecule::Fragmentable do

  let(:test_mol) { "COP(=O)(O)OCNCOCC(OO)C(=O)O" }

  describe 'the :cn rule' do

    it 'cleaves like an ether a secondary NH group if possible' do
      mol = Rubabel["CCNC"]
      frag_sets = mol.fragment(rules: [:cn])
      frag_sets.size.should == 1
      csmiles = frag_sets.first.map(&:csmiles)
      csmiles.should include("C=C")
      csmiles.should include("C[NH3+]")
    end

    it 'will not cleave if not possible' do
    end


  end


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

      it 'lops off carboxy groups' do
        mol = Rubabel["NCC(=O)O"]
        frag_sets = mol.fragment( rules: [:co] )
        frag_sets.size.should == 1
        csmiles = frag_sets.first.map(&:csmiles)
        
        csmiles.should include("[CH2-][NH3+]")
        csmiles.should include("O=C=O")
      end

      it "doesn't remove the group if it is not c3" do
        mol = Rubabel["C=CC(=O)O"]
        fragments = mol.fragment( rules: [:co] )
        fragments.should be_empty
      end

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

      it 'does not cleave esters without position for double bond' do
        mol = Rubabel["NCCC(=O)OC"]
        pieces = mol.fragment( rules: [:co] )
        pieces.should be_empty
      end

      it 'does cleave esters with place for double bond' do
        mol = Rubabel["NCCC(=O)OCC"]
        pieces = mol.fragment( rules: [:co] )
        pieces.size.should == 1  # one set
        the_pair = pieces.first
        csmiles = the_pair.map(&:csmiles)
        csmiles.should include("OC(=O)CC[NH3+]")
        csmiles.should include("C=C")
      end

    end
  end


end
