require 'spec_helper'

require 'pry'
require 'rubabel/molecule'

describe Rubabel::Molecule do
  describe 'breaking a molecule' do
    before(:each) do
      @mol =  Rubabel::Molecule.from_string("NC(=O)CO")
      @n = @mol.find {|a| a.el == :N }
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
      p pieces
      pieces.map(&:csmiles).sort.should == ["N", "OCC=O"]
    end
    describe "Adduct Cases:"

    it 'splits single break with a single adduct' do 
      mol = Rubabel::Molecule.from_string("NC(=O)C[O-].[Na+]")
      n = mol.find {|a| a.el == :N }
      mol.adducts.size.should == 1
      mol.adduct?.should be_true
      mol.delete_bond(n, n.atoms.first)
      pieces = mol.split
      pieces.size.should == 2 # Or 4 as I move the adduct around?
      pieces.map {|a| a.size.should == 2}
    end
    it 'splits single break with multiple adducts' do 
      mol = Rubabel::Molecule.from_string("[K+].NC(=O)C[O-].[Na+]")
      n = mol.find {|a| a.el == :N }
      mol.adducts.size.should == 2
      mol.adduct?.should be_true
      mol.delete_bond(n, n.atoms.first)
      pieces = mol.split
      pieces.size.should == 3 # Or 4 as I move the adduct around?
      pieces.map {|a| a.size.should == 2}
    end
    it 'splits multiple breaks and single adduct' do 
      mol = Rubabel::Molecule.from_string("NC(=O)C[O-].[Na+]")
      num_bonds_before = mol.num_bonds
      num_atoms_before = mol.num_atoms
      reply = mol.split(mol.bonds.first, mol.bonds.last)
      
      reply.should be_a(Array)
      reply.size.should == 3
      mol.num_bonds.should == num_bonds_before
      mol.num_atoms.should == num_atoms_before
      csmiles = reply.map {|a| a.map(&:csmiles) }.flatten
      csmiles.sort.should == %w(N N.[Na+] CC=O CC=O.[Na+] O O.[Na+]).sort

    end
    it 'splits multiple breaks and multiple adducts' do 
      pending "Do we need to handle multiple adducts?"
    end
    it 'handles the product arrays properly' do 
      mol = Rubabel::Molecule.from_string("NC(=O)C[O-].[Na+]")
      reply = mol.split(mol.bonds.first, mol.bonds.last)
      adducts = mol.adducts
      p adducts
      mols = mol.ob.separate.map(&:upcast).delete_if {|a| adducts.include?(a)}
      p mols
      mols2 = mols.map(&:dup)
      adducts.each do |adduct|
        mols2.each {|mol| mol.associate_atom! adduct.atoms.first}
      end
      products = mols.product(mols2)
      p products
      p mol.formula
      p products[1].map(&:formula)
      p products[2].map(&:formula)
    end

    it 'can iterate through fragments' do
      expected = %w(N OCC=O)
      @mol.delete_bond(@n, @n.atoms.first)
      @mol.each_fragment do |frag|
        frag.csmiles.should == expected.shift
      end
    end
  end
end
