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

  describe '::molecule_from_file and ::molecule_from_string' do
    before(:each) do 
      @samples = TESTFILES + "/Samples.sdf"
    end

    it 'return a single molecule (the first one in the file)' do
      mol_f = Rubabel.molecule_from_file(@samples)
      mol_s = Rubabel.molecule_from_string(IO.read(@samples), :sdf)
      [mol_f, mol_s].each {|mol| mol.should be_a(Rubabel::Molecule) }
      mol_f.should == mol_s
    end
  end

  describe 'can deal with .gz files properly' do
    before(:each) do 
      @gz_file = TESTFILES + "/Samples.sdf.gz"
    end

    it 'can get the file format' do
      # non-existant file okay
      Rubabel.format_from_ext("silly.sdf.gz").should == :sdf
      Rubabel.format_from_ext(@gz_file).should == :sdf
      Rubabel.foreach(@gz_file) do |mol|
        mol.should be_a(Rubabel::Molecule)
      end
    end

  end

  describe 'format from extension' do

    it 'determines format from extension' do
      Rubabel.format_from_ext( TESTFILES + "/Samples.sdf" ).should == :sdf
      Rubabel.format_from_ext( TESTFILES + "/Samples.non_existent" ).should be_nil
    end

    it 'determines format from mime-type' do
      Rubabel.format_from_mime( "chemical/x-mdl-sdfile" ).should == :sdf
      Rubabel.format_from_mime( "chemical/wierdness" ).should be_nil
    end

  end

  #describe 'an atom in it' do
  #end
end

