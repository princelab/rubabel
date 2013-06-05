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

    describe 'cod: carbonyl appendage dump' do
      # a primary oxygen or peroxide => C=O appendage dump
      
      specify 'cod: primary alcohol' do
        mol = Rubabel["NCC(O)CC"]
        frags = mol.fragment(rules: [:cod])
        frags.flatten(1).map(&:csmiles).should == ["C[NH3+]", "CCC=O", "C([NH3+])C=O", "CC"]
      end

      specify 'peroxide' do
        mol = Rubabel["NCC(OO)CC"]
        frags = mol.fragment(rules: [:codoo])
        frags.flatten(1).map(&:csmiles).should == ["OC[NH3+]", "CCC=O", "C([NH3+])C=O", "CCO"]
      end

      specify 'cod: carboxylate' do
        mol = Rubabel["CCC(=O)O"]
        pieces = mol.fragment(rules: [:cod])
        pieces.flatten(1).map(&:csmiles).should == ["[CH2-]C", "O=C=O"]
      end

      specify 'cod: carboxylic acid' do
        mol = Rubabel["CCC(=O)O"]
        mol.add_h!(1.5)
        pieces = mol.fragment(rules: [:cod])
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
        #ff.first.formula.should == "C2H5"
        ff.first.formula.should == "C2H5-"
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
        #ff.first.formula.should == "C3H5O"
        ff.first.formula.should == "C3H5O+"
        ff.first.exact_mass.should be_within(1e-6).of(57.034039779909996)
        #ff.last.formula.should == "O"
        ff.last.formula.should == "O-"
      end

      specify 'phosphodiester' do
        mol = Rubabel["CC(COP(=O)([O-])OCCN"]
        frag_set = mol.fragment(rules: [:oxepd])
        ff = frag_set.first
        ff.first.csmiles.should == '[O-]CCC' 
        ff.last.csmiles.should == '[NH3+]CCO[P](=O)=O'
        #ff.first.formula.should == 'C3H7O'
        ff.first.formula.should == 'C3H7O-'
        ff.first.exact_mass.should be_within(1e-6).of(59.049689844)
        #ff.last.formula.should == 'C2H7NO3P'
        ff.last.formula.should == 'C2H7NO3P+'
        ff.last.exact_mass.should be_within(1e-6).of(124.016354719)

        mol = Rubabel["CCCOP(=O)(OCC[N+](C)(C)C)[O-]"]
        frag_set = mol.fragment(rules: [:oxepd, :oxe])
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
    describe 'oxh: oxygen alpha/beta/gamma hydrogen stealing' do
      specify 'primary alcohol giving water loss' do
        mol = Rubabel["CC(O)CCN"]
        frags = mol.fragment(rules: [:oxh])
        ff = frags.first
        ff.first.csmiles.should == 'C=CCC[NH3+]'
        ff.last.csmiles.should == 'O'
        ll = frags.last
        ll.first.csmiles.should == 'CC=CC[NH3+]'
        ll.last.csmiles.should == 'O'
        #ff.first.formula.should == 'C4H10N'
        ff.first.formula.should == 'C4H10N+'
        ff.first.exact_mass.should be_within(1e-6).of(72.0813243255)
      end

      specify 'peroxide carbonyl formation (or peroxide formation [that what we want??])' do
        # do we really see peroxide formation?  Tamil didn't include this in
        # the rules but it follows from the broad way for creating these
        # rules.  Can prohibit peroxide formation in future if necessary...
        mol = Rubabel["CC(OO)CCN"]
        frags = mol.fragment(rules: [:oxh])
        mols = frags.flatten
        mols.map(&:csmiles).should == ["C=CCC[NH3+]", "OO", "CC(=O)CC[NH3+]", "O", "CC=CC[NH3+]", "OO"]
        mols.map(&:formula).should == ["C4H10N+", "H2O2", "C4H10NO+", "H2O", "C4H10N+", "H2O2"]
        #mols.map(&:formula).should == ["C4H10N", "H2O2", "C4H10NO", "H2O", "C4H10N", "H2O2"]
        mols.map(&:exact_mass).zip([72.081324325, 34.005479304, 88.076238945, 18.010564684, 72.081324325, 34.005479304]) do |act, exp|
          act.should be_within(1e-6).of(exp)
        end
      end

      specify 'ether to alcohol, ignoring errors' do
        # this is a good example of a 'disallowed structure' where the
        # formula's do not match up to the original formulas
        mol = Rubabel["CCOCCN"]
        frags = mol.fragment(rules: [:oxh], errors: :ignore)
        mols = frags.flatten
        mols.map(&:csmiles).should == ["C=C", "OCC[NH3+]", "CCO", "C=C[NH2+]"]
      end

      specify 'ether to alcohol, removing errors' do
        mol = Rubabel["CCOCCN"]
        frags = mol.fragment(rules: [:oxh])
        mols = frags.flatten
        mols.map(&:csmiles).should == ["C=C", "OCC[NH3+]"]
      end

      specify 'ester to alcohol' do
        mol = Rubabel["CC(=O)OCCCN"]
        frags = mol.fragment(rules: [:oxh])
        mols = frags.flatten
        mols.map(&:csmiles).should == ["C=C=O", "OCCC[NH3+]", "CC(=O)O", "C=CC[NH3+]"]
        #mols.map(&:formula).should == ["C2H2O", "C3H10NO", "C2H4O2", "C3H8N"]
        mols.map(&:formula).should == ["C2H2O", "C3H10NO+", "C2H4O2", "C3H8N+"]
        mols.map(&:exact_mass).zip([42.010564684, 76.076238945, 60.021129368000004, 58.065674261]) do |act,exp|
          act.should be_within(1e-6).of(exp)
        end
      end

      specify 'phosphodiester (right now needs very low pH and NOT SURE WORKING PROPERLY)' do
        mol = Rubabel["CC(COP(=O)(O)OCCCN"]
        mol.add_h!(1.0)
        frags = mol.fragment(rules: [:oxhpd])
        frags = frags.flatten
        frags.map(&:csmiles).should == ["CCCO", "[NH3+]CCCOP(=O)=O", "CCCOP(=O)=O", "OCCC[NH3+]"]
      end
    end

  end
end



