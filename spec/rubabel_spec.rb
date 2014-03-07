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

  describe "molecule from sdf archive string" do 

    string = "LMFA00000001|  LIPDMAPS04231208332D|| 45 44  0  0  0  0  0  0  0  0999 V2000|   16.4440    9.0927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   17.6562    7.5630    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   16.4440    9.8589    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   15.7757    8.7097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   15.1072    9.0927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   14.4385    8.7097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   13.7699    9.0927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   13.1013    9.4758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   12.4326    9.8589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   11.7642    9.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   11.0955    9.8589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   10.4269    9.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.7583    9.8589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.0896    9.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    8.4210    9.8589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.7524    9.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.0839    9.8589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    6.4152    9.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    5.7466    9.8589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    5.0780    9.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.0896    8.7097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   16.3707    6.0144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   16.3584    5.2483    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   15.7088    6.4081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   15.0340    6.0359    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   14.3717    6.4297    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   13.6970    6.0576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   13.0223    5.6854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   12.3475    5.3132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   11.6852    5.7070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   11.0105    5.3348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   10.3482    5.7285    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.6735    5.3564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.0111    5.7502    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    8.3364    5.3780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    7.6741    5.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    6.9994    5.3997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    6.3370    5.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    5.6623    5.4213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    5.0000    5.8150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|    9.0235    6.5162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   15.7757    7.9437    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   15.1086    7.5585    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|   15.7212    7.1740    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0|   16.3945    7.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0|  2  1  1  0  0  0  0|  1  3  2  0  0  0  0|  1  4  1  0  0  0  0|  4  5  1  0  0  0  0|  5  6  1  0  0  0  0|  6  7  1  0  0  0  0|  7  8  3  0  0  0  0|  8  9  1  0  0  0  0|  9 10  1  0  0  0  0| 10 11  1  0  0  0  0| 11 12  1  0  0  0  0| 12 13  1  0  0  0  0| 13 14  1  0  0  0  0| 14 15  1  0  0  0  0| 15 16  1  0  0  0  0| 16 17  1  0  0  0  0| 17 18  1  0  0  0  0| 18 19  1  0  0  0  0| 19 20  2  0  0  0  0| 14 21  1  0  0  0  0| 22 23  2  0  0  0  0| 22 24  1  0  0  0  0| 24 25  1  0  0  0  0| 25 26  1  0  0  0  0| 26 27  1  0  0  0  0| 27 28  3  0  0  0  0| 28 29  1  0  0  0  0| 29 30  1  0  0  0  0| 30 31  1  0  0  0  0| 31 32  1  0  0  0  0| 32 33  1  0  0  0  0| 33 34  1  0  0  0  0| 34 35  1  0  0  0  0| 35 36  1  0  0  0  0| 36 37  1  0  0  0  0| 37 38  1  0  0  0  0| 38 39  1  0  0  0  0| 39 40  2  0  0  0  0| 34 41  1  0  0  0  0|  2 22  1  0  0  0  0|  4 42  1  0  0  0  0| 42 43  1  0  0  0  0| 24 44  1  0  0  0  0| 44 45  1  0  0  0  0|M  STY  2   1 SUP   2 SUP|M  SAL   1  2  42  43|M  SBL   1  1  41|M  SMT   1 OCH3|M  SAL   2  2  44  45|M  SBL   2  1  43|M  SMT   2 OCH3|M  END"
    mol = Rubabel.molecule_from_archive(string)
    p mol.class
    mol.is_a?(Rubabel::Molecule).should == true
    #mol.should be_a(Rubabel::Molecule)
    mol2 = Rubabel["LMFA00000001", :lmid]
    mol.should == mol2
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

