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

  describe 'manipulating the bond order' do
    subject { 
      mol = Rubabel['CC'] 
      mol[0].get_bond(mol[1])
    }

    specify '#inc! increases the bond order' do
      subject.order.should == 1
      subject.inc! 
      subject.order.should == 2
    end

    specify '#inc!(n) increases the bond order by n' do
      subject.order.should == 1
      subject.inc!(2) 
      subject.order.should == 3
    end

    specify '#inc!(-n) will not go below 0' do
      subject.order.should == 1
      subject.inc!(-1) 
      subject.order.should == 0
      subject.inc!(-1) 
      subject.order.should == 0
    end

    specify '#+(n) will add but not return same ruby object (underlying ob object same)' do
      other = subject + 1
      other.order.should == 2
      subject.order.should == 2
      subject.should_not equal(other)
      #other.object_id.should_not == subject.object_id
    end

    specify '#+(n) will add but not return same ruby object (underlying ob object same)' do
      another_var = subject
      another_var += 1
      another_var.order.should == 2
      subject.order.should == 2
      subject.ob.should equal(another_var.ob)
    end

    specify '#dec!' do
      subject.dec!
      subject.order.should == 0
    end

    specify '#-(n)' do
      var = subject
      var.inc!
      var.order.should == 2
      var -= 2
      var.order.should == 0
      var -= 1
      var.order.should == 0
    end
  end
end
