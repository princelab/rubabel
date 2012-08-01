require 'spec_helper'

require 'rubabel'

describe Rubabel::Molecule::Fragmentable do

  # :peroxy_to_carboxy
  # :oxygen_asymmetric_sp3, :nitrogen_asymmetric_sp3,
  # :internal_phosphoester

  describe 'fragmentation rules' do
    # coenzyme: CC1=CC(=O)C=CC1=O
    # 2-methylcyclohexa-2,5-diene-1,4-dione

    let(:test_mol) { "COP(=O)(O)OCNCOCC(OO)C(=O)O" }

    it 'raises an error for a bad rule' do
      mol = Rubabel["CCNC"]
      expect { mol.fragment(rules: [:wackiness]) }.to raise_error
    end

    describe ':sp3c_nitrogen_double_bond', :pending do

      it 'cleaves like an ether a secondary NH group if possible' do
        mol = Rubabel["CCNC"]
        frag_sets = mol.fragment(rules: [:sp3c_nitrogen_double_bond])
        frag_sets.size.should == 1
        csmiles = frag_sets.first.map(&:csmiles)
        csmiles.should include("C=C")
        csmiles.should include("C[NH3+]")
      end

      it 'will not cleave if not possible' do
        mol = Rubabel["CNC"]
        frag_sets = mol.fragment(rules: [:sp3c_nitrogen_double_bond])
        frag_sets.should be_empty
      end

    end

    describe ':co2_loss', :pending do
      it 'loss of CO2 from carboxy group with charge transfer' do
        mol = Rubabel["NCC(=O)O"]
        frag_sets = mol.fragment( rules: [:co2_loss] )
        frag_sets.size.should == 1
        csmiles = frag_sets.first.map(&:csmiles)

        csmiles.should include("[CH2-][NH3+]")
        csmiles.should include("O=C=O")
      end

      it "doesn't remove CO2 if adjacent is not c3" do
        mol = Rubabel["C=CC(=O)O"]
        fragments = mol.fragment( rules: [:co2_loss] )
        fragments.should be_empty
      end

    end

    describe ':peroxy_to_carboxy' do
      mol = Rubabel["NCCC(OO)CC"]
      fragments = mol.fragment( rules: [:peroxy_to_carboxy] )
    end

    describe ':oxygen_asymmetric_sp3', :pending do
      it 'splits like sp3c_oxygen_double_bond except oxygen takes the electrons' do
        mol = Rubabel["NCCC(O)CC"]
        mol = Rubabel["NCC(O)CC"]
        mol = Rubabel["NC(O)CC"] 
      end
    end

    describe ':sp3c_oxygen_double_bond', :pending do

      describe 'water loss' do
        it 'does h2o loss of alcohol' do
          mol = Rubabel["NCCC(O)CC"]
          fragments = mol.fragment( rules: [:sp3c_oxygen_double_bond] )
          fragments.flatten(1).map(&:csmiles).sort.should == ["CC=CCC[NH3+]", "CCC=CC[NH3+]", "O", "O"]
        end

        it 'h2o loss does not allow bad chemistry' do
          # lone pair and double bond resonance ?
          mol = Rubabel["NCC(O)CC"]
          fragments = mol.fragment( rules: [:sp3c_oxygen_double_bond] )
          fragments.flatten(1).map(&:csmiles).sort.should == ["CC=CC[NH3+]", "O"]

          mol = Rubabel["NC(O)CC"] 
          fragments = mol.fragment( rules: [:sp3c_oxygen_double_bond] )
          fragments.flatten(1).map(&:csmiles).sort.should == []
        end
      end

      describe 'backbone cleavage', :pending do

        it 'does not cleave esters without sp3 carbons available for double bond' do
          mol = Rubabel["NCCC(=O)OC"]
          pieces = mol.fragment( rules: [:sp3c_oxygen_double_bond] )
          pieces.should be_empty
        end

        it 'cleaves esters on far side' do
          mol = Rubabel["NCCC(=O)OCC"]
          pieces = mol.fragment( rules: [:sp3c_oxygen_double_bond] )
          pieces.size.should == 1  # one set
          the_pair = pieces.first
          csmiles = the_pair.map(&:csmiles)
          csmiles.should include("OC(=O)CC[NH3+]")
          csmiles.should include("C=C")
        end

      end
    end

    describe ':alcohol_to_aldehyde', :pending do
      it 'cleaves beside alcohols to generate an aldehyde' do
        mol = Rubabel["NCCC(O)CC"]
        mol.correct_for_ph!
        total_mass = mol.add_h!.mass

        pieces = mol.fragment(rules: [:alcohol_to_aldehyde])
        pieces.size.should == 2
        pieces.map(&:size).should == [2,2]
        pieces.flatten(1).map(&:csmiles).should == ["CC[NH3+]", "CCC=O", "C(C=O)C[NH3+]", "CC"]
        pieces.each do |pair|
          pair.map(&:mass).reduce(:+).should == total_mass
        end
      end
    end

  end
end
