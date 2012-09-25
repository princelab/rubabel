require 'spec_helper'

require 'rubabel/molecule'
require 'rubabel/bond'

describe Rubabel::Bond do
  describe 'cholesterol from sdf' do
    subject { Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' ).bonds.first }

    it 'is a Rubabel::Bond' do
      subject.should be_a(Rubabel::Bond)
    end

    it 'knows what atoms it includes' do
      subject.each_atom do |atom|
        atom.should be_a(Rubabel::Atom)
      end
      subject.atoms.size.should == 2
    end
  end

end
