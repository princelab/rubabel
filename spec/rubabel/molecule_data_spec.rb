require 'spec_helper'

require 'rubabel/molecule_data'
require 'rubabel'

describe Rubabel::MoleculeData do

  before(:each) do
    sdf = TESTFILES + "/Samples.sdf"
    @mol = Rubabel::Molecule.from_file( sdf )
  end

  it 'is initialized with an OpenBabel molecule' do
    md = Rubabel::MoleculeData.new(@mol.ob)
    md.should be_a(Rubabel::MoleculeData)
  end

  describe 'hash-like behavior' do
    subject { Rubabel::MoleculeData.new(@mol.ob) }

    it '#each' do
      enum = subject.each
      enum.should be_a(Enumerator)
      pair = enum.next
      pair.should == ["NAME", "2Ferrocene"]
    end

    it '#to_a' do
      subject.to_a.size.should == 2
    end

    it '#size & #length' do
      subject.size.should == 2
      subject.length.should == 2
    end
    
    it '#[]' do
      subject['NAME'].should == "2Ferrocene"
    end

    it '#[]=' do
      # modify one:
      subject['NAME'] = 'PEPPER'
      subject.size.should == 2
      # create_a new one:
      subject['jtp_special'] = 'sauce'
      subject.size.should == 3
      string = subject.obmol.upcast.write(:sdf)
      string.should =~ /jtp_special/
      string.should =~ /sauce/
      string.should =~ /PEPPER/
    end

    it '#key?' do
      subject.key?('NAME').should be_true
      subject.key?('bananas').should be_false
    end

    it '#keys' do
      subject.keys.should == ["NAME", "OpenBabel Symmetry Classes"]
    end

    it '#values' do
      subject.values.should == ["2Ferrocene", "8 4 9 4 4 4 4 4 4 4 4 4 5 3 6 2 7 1 5 3 6 2 5 3 7 1 5 3 6 2 6 2"]
    end

    it '#delete' do
      key = "OpenBabel Symmetry Classes"
      subject.delete(key).should =~ /8 4 9/
      subject.key?(key).should be_false
      subject.size.should == 1
      subject.delete("nonsense").should be_nil
      subject.delete("nonsense") { 'wow' }.should == 'wow'
    end
  end
end
