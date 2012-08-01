require 'spec_helper'

require 'rubabel/molecule'

describe Rubabel::Molecule do
  describe 'creation' do
    it 'can be made with Rubabel[]' do
      mol = Rubabel["CC(O)O"]
      mol.csmiles.should == "CC(O)O"
    end
  end

  before(:each) do
    @mol = Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' )
  end

  attributes = {
    charge: 0,
    formula: "C27H46O",
    dim: 2,
    spin: 1,
  }.each do |methd, exp|
    it "has #{methd}" do
      @mol.send(methd).should == exp
    end
  end

  specify '#swap! can swap atoms' do
    mol = Rubabel["NCC(=O)O"]
    atoms = mol.atoms
    p mol.atoms
    p mol.bonds
    p mol
    p mol.swap!(atoms[1], atoms[0], atoms[3], atoms[4])
    p mol
    p mol.atoms
    p mol.bonds
  end

  specify '#dup creates an entirely new molecule based on the first' do
    another = @mol.dup
    # this is a deep copy all the way.  Even the atoms are duplicated so that
    # they can be modified in one and do not affect the other at all.
    @mol.atoms.first.charge = 1
    @mol.charge.should_not == another.charge
  end

  specify '#each iterates through each atom in id order' do
    cnt = 0
    @mol.each do |atom|
      atom.id.should == cnt
      cnt += 1
    end
    @mol.atoms.size.should == 33
    @mol.add_h!
    @mol.atoms.size.should == 74
  end

  specify '#hydrogens_added?' do
    @mol.hydrogens_added?.should be_false
    @mol.atoms.size.should == 33
    @mol.add_h!
    @mol.atoms.size.should == 74
    @mol.hydrogens_added?.should be_true
  end

  it 'calculates #ob_sssr (smallest set of smallest rings)' do
    ar = @mol.ob_sssr
    ar.should be_an(Array)
    # in the future should be Rubabel::Ring
    ar.first.should be_a(OpenBabel::OBRing)
  end

  describe 'masses' do
    subject { Rubabel::Molecule.from_string("C(=O)COC(=O)C[NH3+]") }
    it '#mol_wt (or #avg_mass)' do
      subject.mol_wt.should be_within(0.000001).of(118.11121999999999)
    end

    it '#exact_mass' do
      subject.exact_mass.should be_within(0.00000001).of(118.05041812003999)
    end

    it '#mass is the exact mass adjusted for electron gain/loss' do
      subject.mass.should be_within(0.00000001).of(118.04986952003999)
    end
  end

  describe 'getting other descriptors' do
    # can't figure this out yet
  end

  describe 'pH' do

    subject { Rubabel::Molecule.from_string("NCC(=O)OCC(=O)O") }

    it '#correct_for_ph! neutral' do
      subject.correct_for_ph!.to_s.should == '[O-]C(=O)COC(=O)C[NH3+]'
    end

    it '#correct_for_ph!(1.4) [low]' do
      subject.correct_for_ph!(1.4).to_s.should == 'OC(=O)COC(=O)C[NH3+]'
    end

    it '#correct_for_ph!(11.0) [high]' do
      subject.correct_for_ph!(11.0).to_s.should == '[O-]C(=O)COC(=O)CN'
    end

    it '#correct_for_ph!(nil) [gives neutral molecule]' do
      subject.correct_for_ph!.to_s.should == '[O-]C(=O)COC(=O)C[NH3+]'
      subject.correct_for_ph!(nil).to_s.should == "NCC(=O)OCC(=O)O"
    end

    it '#neutral! [can be set neutral again]' do
      subject.correct_for_ph!
      subject.to_s.should == '[O-]C(=O)COC(=O)C[NH3+]'
      subject.h_added?.should == false
      subject.neutral!.to_s.should == "NCC(=O)OCC(=O)O"
      subject.h_added?.should == false
    end

    it '#neutral! [preserves hydrogens added state]' do
      subject.correct_for_ph!
      subject.add_h!
      subject.to_s.should == '[O-]C(=O)COC(=O)C[NH3+]'
      subject.h_added?.should == true
      subject.neutral!.to_s.should == "NCC(=O)OCC(=O)O"
      subject.h_added?.should == true
    end

    it '#add_h!(11.0) [can correct for ph if given a ph]' do
      subject.add_h!(11.0).to_s.should == '[O-]C(=O)COC(=O)CN'
    end

  end

  it 'can calculate a molecular fingerprint (for similarity calcs)' do
    # this just returns the std::vector object at the moment
    fp = @mol.ob_fingerprint
    # this is an array of unsigned ints that really need to be coerced into
    # bits for further usefulness.
    fp.to_a.should == [0, 604110848, 16777216, 0, 2147483648, 4210688, 0, 2097152, 16, 16809984, 0, 0, 1, 37756928, 32, 0, 524296, 1028, 8388612, 131072, 1073741824, 512, 1048584, 16384, 1026, 0, 0, 524288, 0, 2048, 16777248, 0]
    lambda { @mol.ob_fingerprint("WACKY") }.should raise_error
  end

  it 'can calculate the tanimoto similarity' do
    # oo way to call:
    @mol.tanimoto(@mol).should == 1.0
    mol2 = Rubabel::Molecule.from_string("CCC(O)OCCC")
    # class way to call this
    t = Rubabel::Molecule.tanimoto(@mol, mol2)
    # actual: 0.11363636363636363
    t.should be < 0.2
    t.should be > 0.0
  end

  describe '3D' do
    before(:each) do
      # cholesterol
      @mol = Rubabel::Molecule.from_string("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    end

    it 'can be turned into a 3D molecule' do
      # this is only as good as the Builder is good.  For instance, it fails
      # to get all the stereo centers of cholesterol (but it does warn on
      # this, although I don't know how to capture the warnings (can't get
      # with stdout or stderr??))
      @mol.ob.has_3d.should be_false
      @mol.make_3d!
      @mol.ob.has_3d.should be_true 
    end
  end

  describe 'breaking a molecule' do
    before(:each) do
      @mol =  Rubabel::Molecule.from_string("NC(=O)CO")
    end

    it 'num_atoms, atoms and each_atom are sensitive to #add_h!' do
      @mol.num_atoms.should == 5
      @mol.atoms.size.should == 5
      @mol.each_atom.map.to_a.size.should == 5
      @mol.add_h!
      @mol.num_atoms.should == 10
      @mol.atoms.size.should == 10
      @mol.each_atom.map.to_a.size.should == 10
    end

    it 'num_bonds, bonds and each_bond are also sensitive to #add_h!' do
      @mol.num_bonds.should == 4
      @mol.bonds.size.should == 4
      @mol.each_bond.map.to_a.size.should == 4
      @mol.add_h!
      @mol.num_bonds.should == 9
      @mol.bonds.size.should == 9
      @mol.each_bond.map.to_a.size.should == 9
    end

    it 'can be split into multiple molecules [unaffecting self]' do
      num_bonds_before = @mol.num_bonds
      num_atoms_before = @mol.num_atoms

      reply = @mol.split(@mol.bonds.first, @mol.bonds.last)

      reply.should be_a(Array)
      reply.size.should == 3
      @mol.num_bonds.should == num_bonds_before
      @mol.num_atoms.should == num_atoms_before
      csmiles = reply.map(&:csmiles)
      csmiles.sort.should == %w(N CC=O O).sort
    end

  end

  describe 'matching patterns (SMARTS)' do
    before(:each) do
      @mol =  Rubabel::Molecule.from_string("NC(=O)O")
    end

    it 'can match smarts patterns' do
      smarts_pattern = 'C=O'
      ar = @mol.matches(smarts_pattern)
      ar.should be_a(Array)
      ar.size.should == 1
      ar.first.map(&:type).should == %w(Cac O.co2)

      # reverse order of smiles gives reverse order
      ar = @mol.matches(smarts_pattern.reverse)
      ar.first.map(&:type).should == %w(Cac O.co2).reverse
    end

  end
end
