require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/atom'

describe Rubabel::Atom do

  xit 'can be created given an element symbol' do
    hydrogen = Rubabel::Atom[:h]
    hydrogen.el.should == :h

    carbon = Rubabel::Atom[:c]
    carbon.el.should == :c

    chlorine = Rubabel::Atom[:cl, 3]
    chlorine.el.should == :cl
    chlorine.id.should == 3
  end

  describe 'attaching another atom to itself (as part of a molecule)' do

    it 'attaches preformed atom' do
    end

    it 'attaches given an atomic number' do

    end

    it 'attaches given an element symbol' do
    end

  end

  specify 'equality' do
    mol = Rubabel["CCO"]
    oxygen = mol.atoms[2]
    oxygen_from_match = mol.matches("CO").first.last
    (oxygen == oxygen_from_match).should be_true
    (oxygen.equal?(oxygen_from_match)).should be_true
    (oxygen.eql?(oxygen_from_match)).should be_true

    mol2 = Rubabel["CCO"]
    (mol.atoms[0] == mol2.atoms[0]).should_not be_true
    (mol.atoms[0].equal?(mol2.atoms[0])).should_not be_true
    (mol.atoms[0].eql?(mol2.atoms[0])).should_not be_true
  end


  describe 'removing protons and hydrides' do

    def are_identical(obs, exp, atom_index)
      [:formula, :csmiles, :charge, :exact_mass].each do |mol_prop|
        exp_val = exp.send(mol_prop)
        obs_val = obs.send(mol_prop)
        exp_val.is_a?(Float) ? 
          (obs_val.round(5).should == exp_val.round(5)) :
          (obs_val.should == exp_val)
      end
      
      [:formal_charge, :implicit_valence, :isotope, :spin, :valence, :partial_charge].each do |prop|
        obs[atom_index].send(prop).should == exp.atoms[atom_index].send(prop)
      end
    end

    context 'on a molecule with implicit hydrogens' do
      subject { Rubabel["CC"] }

      specify 'removing a hydride' do
        subject[1].remove_a_hydride!
        are_identical( subject, Rubabel["C[CH2+]"].remove_h!, 1)
      end

      specify 'removing a proton' do
        subject[1].remove_a_proton!
        are_identical( subject, Rubabel["C[CH2-]"].remove_h!, 1)
      end
    end

    context 'on a molecule with explicit hydrogens' do
      subject { Rubabel["CC"].add_h! }

      specify 'removing a hydride' do
        subject[1].remove_a_hydride!
        are_identical( subject, Rubabel["C[CH2+]"].add_h!, 1)
      end

      specify 'removing a proton' do
        subject[1].remove_a_proton!
        are_identical( subject, Rubabel["C[CH2-]"].add_h!, 1)
      end

      # the :partial_charge is not matching up on these, not sure why!
      specify 'the atom matches the partial_charge from atom made from smiles'
    end
  end


#  describe 'removing protons and hydrides' do
#describe 'on molecule with implicit hydrogens' do
#subject { Rubabel["CC"] }

      #it 'behaves like an ethane carbanion' do
        ##http://en.wikipedia.org/wiki/Carbanion
        #ethane_carbanion = Rubabel["C[CH2-]"]
        #ethane_carbanion.add_h! if subject.hydrogens_added?
        #[:formula, :csmiles, :charge, :exact_mass].each do |mol_prop|
          #ethane_carbanion.send(mol_prop).should == subject.send(mol_prop)
        #end
        #[:formal_charge, :implicit_valence, :isotope, :partial_charge, :spin, :valence].each do |prop|
          #ethane_carbanion[1].send(prop).should == subject.atoms[altered_atom_idx].send(prop)
        #end
      #end
    #end
  #end

    #xit "removes protons on molecule with implicit hydrogens" do
      #mol = Rubabel["CC"]
      #mol.atoms[1].remove_a_proton!
      #include_examples "is ethane carbanion", mol
    #end

    #xit "removes protons on molecule with explicit H's" do
      #mol = Rubabel["CC"]
      #mol.add_h!
      #mol.atoms[1].remove_a_proton!
      #mol.formula.should == "C2H5-"
      #mol.csmiles.should == '[CH2-]C'
      #mol.exact_mass.round(5).should == 29.03913
      #mol.charge.should == 1
    #end
  #end

  #describe 'properly removing a hydride' do
    #before do
      #@ethane = Rubabel["CC"]
    #end

    #shared_examples "is ethenium" do |mol|
      ## http://en.wikipedia.org/wiki/Ethenium
      ## protonated ethylene
      #ethenium = Rubabel["C[CH2+]"]
      #ethenium.add_h! if mol.hydrogens_added?
      #[:formula, :csmiles, :charge, :exact_mass].each do |mol_prop|
        #ethenium.send(mol_prop).should == mol.send(mol_prop)
      #end
      #[:formal_charge, :implicit_valence, :isotope, :partial_charge, :spin, :valence].each do |prop|
        #ethenium.send(prop).should == mol.send(prop)
      #end
    #end

    #it "removes hydrides on molecule with implicit hydrogens" do
      #@ethane.atoms[1].remove_a_hydride!
      #include_examples "is ethenium", @ethane
    #end

    #xit "removes hydrides on molecule with explicit hydrogens" do
      #@ethane = Rubabel["CC"]
      #@ethane.add_h!
      #@ethane.atoms[1].remove_a_hydride!
      #@ethane.formula.should == "C2H5-"
      #@ethane.csmiles.should == 'C[CH2-]'
      #@ethane.exact_mass.round(5).should == 29.03913
      #@ethane.charge.should == -1
    #end
  #end



  it 'can find atom identities with simple questions' do
    mol = Rubabel["NCC(O)CC(=O)"]
    (c_exp, o_exp) = mol.matches("C=O").first
    mol.find(&:carbonyl_carbon?).id.should == c_exp.id
    mol.find(&:carbonyl_oxygen?).id.should == o_exp.id
  end

  describe 'working with a complex molecule' do

    before do
      @mol = Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' )
      @atom = @mol.atoms.first
      @mol_h =  Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' )
      @mol_h.add_h!
      @atom_with_h = @mol_h.atoms.first
    end

    attributes = {
      charge: 0,
      id: 0,
      spin: 0,
    }.each do |methd, exp|
      it "has #{methd}" do
        @atom.send(methd).should == exp
      end
    end

    specify '#<<(atom) adds an atom and returns the added atom' do
      mol = Rubabel["C"]
      reply = (mol[0] << :n << :o)
      reply.should be_a(Rubabel::Atom)
      mol.csmiles.should == 'CNO'
    end

    it '#mol retrieves the parent molecule' do
      @atom.mol.should == @mol

      ## no parent molecule 
      #h = Rubabel::Atom[:h]
      #h.mol.should be_nil
    end

    it 'can get the bonds' do
      @atom.each_bond do |bond|
        bond.should be_a(Rubabel::Bond)
      end
      @atom.bonds.size.should == 4
    end

    it 'can add a bond'

    it 'can get the neighboring atoms' do
      @atom.id.should == 0
      @atom.atomic_num.should == 6
      @atom.type.should == 'C3'
      @atom.each_atom do |nbr_atom|
        nbr_atom.should be_a(Rubabel::Atom)
        nbr_atom.ob.equal?(@atom.ob).should be_false
      end
      @atom.atoms.size.should == 4
    end

    it '#get_bond can retrieve a particular bond based on the atom' do
      other = @mol.atoms[3]
      bond = @atom.get_bond(other)
      bond.atoms.map(&:id).should == [0,3]
      other = @mol.atoms[15]  # these are not connected
      @atom.get_bond(other).should be_nil
    end

    it '#coords gets the coordinates' do
      @atom.coords
    end
  end
end
