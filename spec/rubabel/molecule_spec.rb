require 'spec_helper'

require 'rubabel/molecule'

describe Rubabel::Molecule do
  describe 'a very simple smiles molecule' do
    subject { Rubabel::Molecule.from_string("CCC(=O)O") }

    it 'has basic characteristics' do
      subject.charge.should == 0
      subject.formula.should == "C3H6O2"
      subject.dim.should == 0
      p subject.num_atoms
      puts "HAER"
      p subject.conformers
      puts 'he'
    end
  end

  describe 'a more complex molecule' do
    subject { Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' ) }

    it 'has basic characteristics' do
      subject.charge.should == 0
      subject.formula.should == "C27H46O"
      subject.dim.should == 2
    end
  end


  
end
