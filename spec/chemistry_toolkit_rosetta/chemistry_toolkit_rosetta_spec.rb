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
    #http://ctr.wikia.com/wiki/Heavy_atom_counts_from_an_SD_file
    output = capture_stdout do
      require 'rubabel'
      Rubabel.foreach("benzodiazepine.sdf.gz") {|mol| puts mol.ob.num_hvy_atoms }
    end
    output.should == IO.read( @keydir + "/" + "benzodiazepine_heavy_atom_counts.output.10.txt" )
  end

  it 'gets Ring counts in a SMILES file' do
    # http://ctr.wikia.com/wiki/Ring_counts_in_a_SMILES_file
    output = capture_stdout do
      require 'rubabel'
      Rubabel.foreach("benzodiazepine.sdf.gz") {|mol| puts mol.ob_sssr.size }
    end
    output.should == IO.read( @keydir + '/' + "benzodiazepine_ring_counts.output.10.txt" )
  end

  it 'Converts a SMILES string to canonical SMILES' do
    # http://ctr.wikia.com/wiki/Convert_a_SMILES_string_to_canonical_SMILES
    require 'rubabel'
    smiles = %w{CN2C(=O)N(C)C(=O)C1=C2N=CN1C CN1C=NC2=C1C(=O)N(C)C(=O)N2C}
    cans = smiles.map {|smile| Rubabel::Molecule.from_string(smile) }
    fail unless cans.reduce(:==)
  end

  it 'Works with SD tag data' do
    # http://ctr.wikia.com/wiki/Working_with_SD_tag_data
    require 'rubabel'
    File.open("RULE5.sdf", 'w') do |out|
      Rubabel.foreach("benzodiazepine.sdf.gz") do |mol|
        mol.data["RULE5"] = 
          if mol.data.key?('PUBCHEM_XLOGP3')
            num_true = [
              mol.data["PUBCHEM_CACTVS_HBOND_DONOR"].to_i <= 5,
              mol.data["PUBCHEM_CACTVS_HBOND_ACCEPTOR"].to_i <= 10,
              mol.data["PUBCHEM_MOLECULAR_WEIGHT"].to_f <= 500,
              mol.data["PUBCHEM_XLOGP3"].to_f <= 5 
            ].count(true)
            if num_true >= 3 then "1" 
            else "0"
            end
          else
            "no logP"
          end
        out.print mol.write(:sdf)
      end
    end
    # end wiki code
    (lines_got, lines_exp) = ["RULE5.sdf", @keydir + "/rule5.10.sdf"].map do |file|
      lines = IO.readlines(file).reject {|line| line =~ /OpenBabel/ }
    end
    lines_got.should == lines_exp
  end
end
