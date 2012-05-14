require 'spec_helper'

# http://ctr.wikia.com/wiki/Chemistry_Toolkit_Rosetta_Wiki


# these are written as if the user is starting the project from scratch in the
# style of the wiki.  We put the test files in the same folder as our specs so
# that the code can be copied directly into the wiki

# the code that should be copied is within the capture_stdout blocks

describe 'Chemistry Toolkit Rosetta Wiki' do

  before(:each) do
    @orig_dir = Dir.pwd
    @wiki_spec_dir = File.dirname(__FILE__)
    Dir.chdir(@wiki_spec_dir)
    @keydir =  @wiki_spec_dir + "/key"
  end

  after(:each) do 
    Dir.chdir(@orig_dir)
  end

  it 'gets Heavy atom counts from an SD file' do
    output = capture_stdout do
      require 'rubabel'
      Rubabel.foreach("benzodiazepine.sdf.gz") {|mol| puts mol.ob.num_hvy_atoms }
    end
    output.should == IO.read( @keydir + "/" + "benzodiazepine_heavy_atom_counts.output.10.txt" )
  end

end
