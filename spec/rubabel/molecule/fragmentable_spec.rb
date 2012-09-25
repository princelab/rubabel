require 'spec_helper'

require 'rubabel'

$VERBOSE = nil

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

    describe 'cad_o: carbonyl appendage dump ' do
      # a primary oxygen or peroxide => C=O appendage dump
      
      describe 'cad_o: primary alcohol' do
        mol = Rubabel["NCC(O)CC"]
        frags = mol.fragment(rules: [:cad_o])
        frags.flatten(1).map(&:csmiles).should == ["C[NH3+]", "CCC=O", "C([NH3+])C=O", "CC"]
      end

      describe 'peroxide' do
        mol = Rubabel["NCC(OO)CC"]
        frags = mol.fragment(rules: [:cad_oo])
        frags.flatten(1).each_with_index do |f,i|
          f.write("mol#{i}.svg")
        end
        frags.flatten(1).map(&:csmiles).should == ["OC[NH3+]", "CCC=O", "C([NH3+])C=O", "CCO"]
      end

      describe 'cad_o: carboxylate' do
        mol = Rubabel["CCC(=O)O"]
        pieces = mol.fragment(rules: [:cad_o])
        pieces.flatten(1).map(&:csmiles).should == ["[CH2-]C", "O=C=O"]
      end

      describe 'cad_o: carboxylic acid' do
        mol = Rubabel["CCC(=O)O"]
        mol.add_h!(1.5)
        pieces = mol.fragment(rules: [:cad_o])
        pieces.flatten(1).map(&:csmiles).should == ["CC", "O=C=O"]
      end
    end

    describe 'oxygen bond stealing' do
      # oxygen just steals the electron pair it is attached to.  This
      # typically results in a negatively charged oxygen and a positively
      # charged carbo-cation.
      describe 'ether to ions' do
      end

      describe 'ester to ions' do
      end

      describe 'carboxyl group' do
      end

      describe 'phosphodiester' do
      end
    end

    # this is really a subset of oxygen bond stealing: if the negatively
    # charged oxygen can rip off a nearby proton, it will.
    describe 'oxygen alpha/beta/gamma hydrogen stealing' do
      describe 'primary alcohol giving water loss' do
      end

      describe 'peroxide carbonyl formation' do
      end

      describe 'ether to alcohol' do
      end

      describe 'ester to alcohol' do
      end

      describe 'phosphodiester' do
      end
    end

  end
end




    #describe ':sp3c_nitrogen_double_bond' do

      #it 'cleaves like an ether a secondary NH group if possible' do
        #mol = Rubabel["CCNC"]
        #frag_sets = mol.fragment(rules: [:sp3c_nitrogen_double_bond])
        #frag_sets.size.should == 1
        #csmiles = frag_sets.first.map(&:csmiles)
        #csmiles.should include("C=C")
        #csmiles.should include("C[NH3+]")
      #end

      #it 'will not cleave if not possible' do
        #mol = Rubabel["CNC"]
        #frag_sets = mol.fragment(rules: [:sp3c_nitrogen_double_bond])
        #frag_sets.should be_empty
      #end

    #end

    #describe ':co2_loss' do
      #it 'loss of CO2 from carboxy group with charge transfer' do
        #mol = Rubabel["NCC(=O)O"]
        #frag_sets = mol.fragment( rules: [:co2_loss] )
        #frag_sets.size.should == 1
        #csmiles = frag_sets.first.map(&:csmiles)

        #csmiles.should include("[CH2-][NH3+]")
        #csmiles.should include("O=C=O")
      #end

      #it "doesn't remove CO2 if adjacent is not c3" do
        #mol = Rubabel["C=CC(=O)O"]
        #fragments = mol.fragment( rules: [:co2_loss] )
        #fragments.should be_empty
      #end

    #end

    #describe ':peroxy_to_carboxy' do
      #it 'works' do
        #mol = Rubabel["NCCC(OO)CC"]
        #frag_sets = mol.fragment( rules: [:peroxy_to_carboxy] )
        #frag_sets.size.should == 2 
        #frag_sets.flatten(1).map(&:csmiles).sort.should == ["CC", "CCC(=O)O", "CC[NH3+]", "OC(=O)CC[NH3+]"]
      #end
    #end

    #describe ':sp3c_oxygen_asymmetric_far_sp3', :pending do
      #it 'splits like sp3c_oxygen_double_bond except oxygen takes the electrons' do
        #$VERBOSE = 3
        #mol = Rubabel["NCCCOCC"]
        #frag_sets = mol.fragment( rules: [:sp3c_oxygen_asymmetric_far_sp3] )
        #$VERBOSE = nil
        #frag_sets.size.should == 2
        ##mol = Rubabel["NCCOCC"]
        ##p mol.fragment( rules: [:sp3c_oxygen_asymmetric_far_sp3] )
        ##mol = Rubabel["NCOC"] 
        ##p mol.fragment( rules: [:sp3c_oxygen_asymmetric_far_sp3] )
      #end
    #end

    #describe ':sp3c_oxygen_double_bond_water_loss' do

      #it 'does h2o loss of alcohol' do
        #mol = Rubabel["NCCC(O)CC"]
        #fragments = mol.fragment( rules: [:sp3c_oxygen_double_bond_water_loss] )
        #fragments.flatten(1).map(&:csmiles).sort.should == ["CC=CCC[NH3+]", "CCC=CC[NH3+]", "O", "O"]
      #end

      #it 'h2o loss does not allow bad chemistry' do
        ## lone pair and double bond resonance ?
        #mol = Rubabel["NCC(O)CC"]
        #fragments = mol.fragment( rules: [:sp3c_oxygen_double_bond_water_loss] )
        #fragments.flatten(1).map(&:csmiles).sort.should == ["CC=CC[NH3+]", "O"]

        #mol = Rubabel["NC(O)CC"] 
        #fragments = mol.fragment( rules: [:sp3c_oxygen_double_bond_water_loss] )
        #fragments.flatten(1).map(&:csmiles).sort.should == []
      #end
    #end

    #describe 'sp3c_oxygen_double_bond_far_side_sp2' do

      #it 'does not cleave esters without sp3 carbons available for double bond' do
        #mol = Rubabel["NCCC(=O)OC"]
        #pieces = mol.fragment( rules: [:sp3c_oxygen_double_bond_far_side_sp2] )
        #pieces.should be_empty
      #end

      #it 'cleaves esters on far side of singly bonded oxygen' do
        #mol = Rubabel["NCCC(=O)OCC"]
        #pieces = mol.fragment( rules: [:sp3c_oxygen_double_bond_far_side_sp2] )
        #pieces.size.should == 1  # one set
        #the_pair = pieces.first
        #csmiles = the_pair.map(&:csmiles)
        #csmiles.should include("OC(=O)CC[NH3+]")
        #csmiles.should include("C=C")
      #end

    #end

    #describe ':alcohol_to_aldehyde' do
      #it 'cleaves beside alcohols to generate an aldehyde' do
        #mol = Rubabel["NCCC(O)CC"]
        #mol.correct_for_ph!
        #total_mass = mol.add_h!.mass

        #pieces = mol.fragment(rules: [:alcohol_to_aldehyde])
        #pieces.size.should == 2
        #pieces.map(&:size).should == [2,2]
        #pieces.flatten(1).map(&:csmiles).should == ["CC[NH3+]", "CCC=O", "C(C=O)C[NH3+]", "CC"]
        #pieces.each do |pair|
          #pair.map(&:mass).reduce(:+).should == total_mass
        #end
      #end
    #end
