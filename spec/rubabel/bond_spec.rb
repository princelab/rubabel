require 'spec_helper'

require 'rubabel/bond'

describe Rubabel::Bond do
  subject { Rubabel::Molecule.from_file( TESTFILES + '/cholesterol.sdf' ).bonds.first }
end
