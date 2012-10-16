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

    #let(:test_mol) { "COP(=O)(O)OCNCOCC(OO)C(=O)O" }

    it 'raises an error for a bad rule' do
      mol = Rubabel["CCNC"]
      expect { mol.fragment(rules: [:wackiness]) }.to raise_error
    end

    describe 'cad_o: carbonyl appendage dump ' do
      # a primary oxygen or peroxide => C=O appendage dump
      
      specify 'cad_o: primary alcohol' do
        mol = Rubabel["NCC(O)CC"]
        frags = mol.fragment(rules: [:cad_o])
        frags.flatten(1).map(&:csmiles).should == ["C[NH3+]", "CCC=O", "C([NH3+])C=O", "CC"]
      end

      specify 'peroxide' do
        mol = Rubabel["NCC(OO)CC"]
        frags = mol.fragment(rules: [:cad_oo])
        frags.flatten(1).map(&:csmiles).should == ["OC[NH3+]", "CCC=O", "C([NH3+])C=O", "CCO"]
      end

      specify 'cad_o: carboxylate' do
        mol = Rubabel["CCC(=O)O"]
        pieces = mol.fragment(rules: [:cad_o])
        pieces.flatten(1).map(&:csmiles).should == ["[CH2-]C", "O=C=O"]
      end

      specify 'cad_o: carboxylic acid' do
        mol = Rubabel["CCC(=O)O"]
        mol.add_h!(1.5)
        pieces = mol.fragment(rules: [:cad_o])
        pieces.flatten(1).map(&:csmiles).should == ["CC", "O=C=O"]
      end
    end

    describe 'oxe: oxygen electron stealing' do
      # oxygen just steals the electron pair it is attached to.  This
      # typically results in a negatively charged oxygen and a positively
      # charged carbo-cation.
      specify 'ether to ions (same for esters)' do
        mol = Rubabel["CCOCCN"]
        frag_set = mol.fragment(rules: [:oxe])
        frags = frag_set.first
        frags.first.csmiles.should == "C[CH2+]"
        frags.last.csmiles.should == '[O-]CC[NH3+]'
        frags.first.formula.should == 'C2H5'
        frags.last.formula.should == 'C2H7NO'
        frags.first.exact_mass.should be_within(1e-6).of(29.03912516)
        frags.last.exact_mass.should be_within(1e-6).of(61.052763849)

        mol = Rubabel["CCOC(=O)CCN"]
        frag_set = mol.fragment(rules: [:oxe])
        ff = frag_set.first
        ff.first.csmiles.should == 'C[CH2+]'
        ff.last.csmiles.should == '[O-]C(=O)CC[NH3+]'
        ff.first.formula.should == "C2H5"
        ff.last.formula.should == "C3H7NO2"
        ff.first.exact_mass.should be_within(1e-6).of(29.03912516035)
        ff.last.exact_mass.should be_within(1e-6).of(89.04767846841)
      end

      specify 'carboxyl group' do
        mol = Rubabel["CCC(=O)O"]
        frag_set = mol.fragment(rules: [:oxe])
        ff = frag_set.first
        ff.first.csmiles.should == 'CC[C+]=O'
        ff.last.csmiles.should == '[O-]'
        ff.first.formula.should == "C3H5O"
        ff.first.exact_mass.should be_within(1e-6).of(57.034039779909996)
        ff.last.formula.should == "O"
      end

      specify 'phosphodiester' do
        mol = Rubabel["CC(COP(=O)([O-])OCCN"]
        frag_set = mol.fragment(rules: [:oxe_phde])
        ff = frag_set.first
        ff.first.csmiles.should == '[O-]CCC' 
        ff.last.csmiles.should == '[NH3+]CCO[P](=O)=O'
        ff.first.formula.should == 'C3H7O'
        ff.first.exact_mass.should be_within(1e-6).of(59.049689844)
        ff.last.formula.should == 'C2H7NO3P'
        ff.last.exact_mass.should be_within(1e-6).of(124.016354719)

        mol = Rubabel["CCCOP(=O)(OCC[N+](C)(C)C)[O-]"]
        frag_set = mol.fragment(rules: [:oxe_phde, :oxe])
        # some of these don't like right on first inspection, but that is
        # because we 'converted dative bonds' meaning + and - next to each
        # other are allowed to cancel one another out!
        frag_set.size.should == 4
        mols = frag_set.flatten
        mols.map(&:csmiles).should ==  ["CC[CH2+]", "[O-]P(=O)(OCC[N+](C)(C)C)[O-]", "CCCOP(=O)([O-])[O-]", "[CH2+]C[N+](C)(C)C", "[O-]CCC", "O=[P](=O)OCC[N+](C)(C)C", "CCCO[P](=O)=O", "[O-]CC[N+](C)(C)C"]
        mols.map(&:formula).should == ["C3H7", "C5H13NO4P", "C3H7O4P", "C5H13N", "C3H7O", "C5H13NO3P", "C3H7O3P", "C5H13NO"]
        mols.map(&:exact_mass).zip([43.05477522449, 182.05821952995, 138.00819533273, 87.10479942171, 59.04968984405, 166.06330491039, 122.01328071317, 103.09971404127]) do |act, exp|
          act.should be_within(1e-6).of(exp)
        end

      end
    end

    # this is really a subset of oxygen bond stealing: if the negatively
    # charged oxygen can rip off a nearby proton, it will.
    describe 'ohs: oxygen alpha/beta/gamma hydrogen stealing' do
      specify 'primary alcohol giving water loss' do
        mol = Rubabel["CC(O)CCN"]
        frags = mol.fragment(rules: [:ohs])
        ff = frags.first
        ff.first.csmiles.should == 'C=CCC[NH3+]'
        ff.last.csmiles.should == 'O'
        ll = frags.last
        ll.first.csmiles.should == 'CC=CC[NH3+]'
        ll.last.csmiles.should == 'O'
        ff.first.formula.should == 'C4H10N'
        ff.first.exact_mass.should be_within(1e-6).of(72.0813243255)
      end

      specify 'peroxide carbonyl formation carbonyl formation or peroxide formation' do
        # do we really get peroxide formation?  Tamil doesn't include this but
        # it follows from the broad way for creating these rules.  Can
        # prohibit peroxide formation in future if necessary...
        mol = Rubabel["CC(OO)CCN"]
        frags = mol.fragment(rules: [:ohs])
        p frags
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
