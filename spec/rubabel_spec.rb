require 'spec_helper'

require 'rubabel'

describe Rubabel do
  describe '::foreach' do
    before(:each) do 
      #lines = %w(7-oxocholesterol-d7.mol cholesterol.mol).flat_map {|fn| IO.readlines(TESTFILES + "/" + fn) + ["\n"] }
      #@two_mols = TESTFILES + '/two.tmp.mol'
      #IO.write( @two_mols, lines.join )
      @two_mols = TESTFILES + "/two.sdf"
    end

    #after(:each) do 
    #  File.unlink @two_mols
    #end

    it 'iterates over every molecule in a file' do
      Rubabel.foreach(@two_mols) do |mol|
        mol.should be_a(Rubabel::Molecule)
      end
      Rubabel.foreach(@two_mols).map.to_a.size.should == 2
    end

    it 'returns an enumerator w/o a block for map or select' do
      Rubabel.foreach(@two_mols).map.should be_a(Enumerator)
    end
  end

  describe '::read_file and ::read_string' do
    before(:each) do 
      @two_mols = TESTFILES + "/two.sdf"
    end

    it 'return a single molecule, the first one' do
      mol = Rubabel.read_file(@two_mols)

      mol = Rubabel.read_string(IO.read(@two_mols))
    end
  end

  #describe 'an atom in it' do
  #end
end

