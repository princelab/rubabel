require 'spec_helper'

require 'rubabel'

describe Rubabel do
  
  def is_format_hash(hash, min_cnt)
    hash = Rubabel.in_formats
    hash.should be_a(Hash)
    # just a ballpark estimate
    hash.size.should be >= min_cnt
    hash.each.next.first.should be_a(Symbol)
  end

  describe '::in_formats' do
    
    it 'gives a hash keyed by symbol' do
      inf = Rubabel.in_formats
      is_format_hash(inf, 120)
    end
  end

  describe '::out_formats' do
    it 'gives a hash keyed by symbol' do
      outf = Rubabel.out_formats
      is_format_hash(outf, 110)
    end
  end

  describe '::foreach' do
    before(:each) do 
      @samples = TESTFILES + "/Samples.sdf"
    end

    it 'iterates over every molecule in a file' do
      cnt = 0
      Rubabel.foreach(@samples) do |mol|
        mol.should be_a(Rubabel::Molecule)
        mol.atoms.size.should be > 0
        cnt += 1
      end
      cnt.should == 41
      Rubabel.foreach(@samples).map.to_a.size.should == 41
    end

    it 'returns an enumerator w/o a block for map or select' do
      Rubabel.foreach(@samples).map.should be_a(Enumerator)
    end
  end

  describe '::read_file and ::read_string' do
    before(:each) do 
      @samples = TESTFILES + "/Samples.sdf"
    end

    it 'return a single molecule (the first one in the file)' do
      mol_f = Rubabel.read_file(@samples)
      mol_s = Rubabel.read_string(IO.read(@samples), :sdf)
      [mol_f, mol_s].each {|mol| mol.should be_a(Rubabel::Molecule) }
      mol_f.should == mol_s
    end
  end

  #describe 'an atom in it' do
  #end
end

