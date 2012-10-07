require 'spec_helper'

require 'rubabel/molecule'

describe Rubabel::Molecule do
  describe 'creation' do
    it 'can be made with Rubabel[]' do
      mol = Rubabel["CC(O)O"]
      mol.csmiles.should == "CC(O)O"
    end

    it 'can be made with an existing obmol object' do
      mol = Rubabel["CC(O)O"]
      mol_dup = Rubabel::Molecule.new(mol.ob)
    end

    it 'can be made fresh' do
      mol = Rubabel::Molecule.new
      mol.atoms.size == 0
      mol.title.should == ''
    end
  end

  #xit 'can add a hydrogen to the formula' do
    #mol = Rubabel["CCC"]
    #p mol.formula
    #mol.add_hydrogen_to_formula!
    #p mol.formula
  #end

  describe 'png output' do
    it 'creates a png image (corresponds to the svg)' do
      mol = Rubabel["NCC(=O)O"]
      #mol.write("mol1.svg", size: '300x200') 
      #mol.write("mol1.png", size: '300x280')
      mol.write("mol1.svg")
      mol.write("mol1.png")
      %w(mol1.svg mol1.png).each do |file|
        File.exist?(file).should be_true
        File.unlink(file)
      end
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

  specify 'eql? and equal? mean the objects modify the same underlying openbabel molecule data' do
    mol = Rubabel["C"]
    eq_mol = mol.atoms.first.mol
    mol.equal?(eq_mol).should be_true
    another = Rubabel["C"]
    mol.equal?(another).should be_false
  end

  specify '== means the canonical smiles strings (:csmiles) are equal' do
    mol1 = Rubabel["CCO"]
    mol2 = Rubabel["OCC"]
    (mol1 == mol2).should be_true
    mol2[0].charge += 1
    (mol1 == mol2).should be_false
    mol3 = Rubabel["CCCO"]
    (mol1 == mol3).should be_false
  end

  specify '#[] retrieves atom by index' do
    mol = Rubabel["NCO"]
    mol[0].el.should == :n
    mol[1].el.should == :c
    mol[2].el.should == :o
    mol[-1].el.should == :o
    mol[-2].el.should == :c
    ar = mol[1..-1]
    ar.first.el.should == :c
    ar.last.el.should == :o
    ar.size.should == 2
  end

  describe 'adding an atom' do
    it 'can be added but not attached' do
      mol = Rubabel["CCO"]
      atom = mol.add_atom!(:n)
      atom.el.should == :n
      atom = mol.add_atom!(8)
      atom.el.should == :o
      mol.csmiles.should == "CCO.N.O"
    end

    it "can be added and attached by el symbol or atomic number" do 
      mol = Rubabel["CCO"]
      first_carbon = mol[0]
      mol.add_atom!(:n, first_carbon)
      mol.csmiles.should == "NCCO"

      mol.add_atom!(16, first_carbon)
      mol.csmiles.should == "NC(CO)S"
    end

  end
  
  specify '#<< adds atom to the last atom and returns the mol' do
    mol = Rubabel::Molecule.new
    # by element symbol or atomic number
    reply = (mol << :n << :c)
    reply.should be_a(Rubabel::Molecule)
    reply.csmiles.should == 'CN'
  end

  specify '#dup duplicates the molecule' do
    mol = Rubabel["CCO"]
    dup_mol = mol.dup
    mol[0].ob.set_atomic_num(9)
    mol.csmiles.should == "OCF"
    dup_mol.csmiles.should == "CCO"
  end

  specify '#add_atom! with a Rubabel::Atom'

  describe '#add_bond! adds a bond (and updates atoms)' do
    specify 'given two atoms' do
      mol = Rubabel["CCO"]
      atom = mol.add_atom!(0)
      mol.add_bond!(mol[1], atom)
      mol.csmiles.should == '*C(O)C'
    end
  end

  specify '#atom(id) retrieves atom by id num' do
    mol = Rubabel["CCO"]
    o = mol.find {|a| a.el == :o }
    mol.atom(o.id).id.should == o.id
  end

  specify '#swap! can swap atoms around' do
    mol = Rubabel["NCC(=O)O"]
    carboxy_carbon = mol.atoms.find {|atom| atom.type == 'Cac' }
    single_bonded_oxygen = mol.atoms.find {|atom| atom.type == 'O.co2' && atom.get_bond(carboxy_carbon).bond_order == 1 }
    nit_carbon = mol.atoms.find {|atom| atom.atoms.any? {|nbr| nbr.type == 'N3' } }
    nitrogen = mol.atoms.find {|atom| atom.el == :n }
    swapped = mol.swap!(nit_carbon, nitrogen, carboxy_carbon, single_bonded_oxygen)
    swapped.should == mol
    swapped.csmiles.should == 'NC(=O)CO'
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

  describe 'making carbo-cations: spin_multiplicity and charges' do
    # http://openbabel.org/docs/2.3.1/Features/Radicals.html
    subject { mol = Rubabel["CC"] }
    it 'can be turned into a carbocation' do
      mol = subject
      c = mol[0]
      c.ob.set_spin_multiplicity 2
      c.charge += 1
      mol.csmiles.should == "C[CH2+]"
    end
  end

  describe 'masses' do
    subject { Rubabel["C(=O)COC(=O)C[NH3+]"] }
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

    subject { Rubabel["NCC(=O)OCC(=O)O"] }

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
      @n = @mol.find {|a| a.el == :n }
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

    it 'can split fragments (akin to separate)' do
      @mol.delete_bond(@n, @n.atoms.first)
      pieces = @mol.split
      pieces.map(&:csmiles).sort.should == ["N", "OCC=O"]
    end

    it 'can iterate through fragments' do
      expected = %w(N OCC=O)
      @mol.delete_bond(@n, @n.atoms.first)
      @mol.each_fragment do |frag|
        frag.csmiles.should == expected.shift
      end
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

  describe 'using output format options' do

  end
end
