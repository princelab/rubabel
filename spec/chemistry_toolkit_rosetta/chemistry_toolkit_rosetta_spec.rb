require 'spec_helper'

# http://ctr.wikia.com/wiki/Chemistry_Toolkit_Rosetta_Wiki


# these are written as if the user is starting the project from scratch in the
# style of the wiki.  We put the test files in the same folder as our specs so
# that the code can be copied directly into the wiki

# the code that should be copied is within the capture_stdout blocks

module Kernel
  def wiki_code(&block)
    block.call
  end

  alias_method :wiki_code_capture_stdout, :capture_stdout
end

def unlink(file)
  File.delete(file) if File.exist?(file)
end

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

  specify 'Heavy atom counts from an SD file' do
    #http://ctr.wikia.com/wiki/Heavy_atom_counts_from_an_SD_file
    output = wiki_code_capture_stdout do
      require 'rubabel'
      Rubabel.foreach("benzodiazepine.sdf.gz") {|mol| puts mol.ob.num_hvy_atoms }
    end
    output.should == IO.read( @keydir + "/" + "benzodiazepine_heavy_atom_counts.output.10.txt" )
  end

  specify 'Ring counts in a SMILES file' do
    # http://ctr.wikia.com/wiki/Ring_counts_in_a_SMILES_file
    output = wiki_code_capture_stdout do
      require 'rubabel'
      Rubabel.foreach("benzodiazepine.sdf.gz") {|mol| puts mol.ob_sssr.size }
    end
    output.should == IO.read( @keydir + '/' + "benzodiazepine_ring_counts.output.10.txt" )
  end

  specify 'Convert a SMILES string to canonical SMILES' do
    # http://ctr.wikia.com/wiki/Convert_a_SMILES_string_to_canonical_SMILES
    wiki_code do
      require 'rubabel'
      smiles = %w{CN2C(=O)N(C)C(=O)C1=C2N=CN1C CN1C=NC2=C1C(=O)N(C)C(=O)N2C}
      cans = smiles.map {|smile| Rubabel[smile] }
      fail unless cans.reduce(:==)
    end
  end

  specify 'Working with SD tag data' do
    # http://ctr.wikia.com/wiki/Working_with_SD_tag_data
    wiki_code do
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
      # NOTE for wiki: This updates the #2 line of each molfile to mention "OpenBabel" and the current timestamp, as it ought to.
    end
    (lines_got, lines_exp) = ["RULE5.sdf", @keydir + "/rule5.10.sdf"].map do |file|
      lines = IO.readlines(file).reject {|line| line =~ /OpenBabel/ }
    end
    lines_got[0,10].should == lines_exp[0,10]
    unlink "RULE5.sdf"
  end

  specify 'Detect and report SMILES and SDF parsing errors' do
    # http://ctr.wikia.com/wiki/Detect_and_report_SMILES_and_SDF_parsing_errors
    wiki_code do
      require 'rubabel'
      File.open("log.txt", 'w') do |out|
        %w(Q C1C).each do |smile|
          Rubabel[smile] rescue out.puts "bad smiles #{smile}"
        end
      end
      # note: error catching can only be performed with the from_string
      # method at the moment.
    end
    IO.read("log.txt").should == "bad smiles Q\nbad smiles C1C\n"
    puts "^^^^^ (ignore the above warning, part of spec) ^^^^^"
    # TODO: implement error catching in file reading
    # TODO: muffle the warning that Open Babel spits out on unmatched ring
    # bonds.  Tried to capture stdout and stderr to no avail.
    unlink("log.txt")
  end

  specify 'Report how many SD file records are within a certain molecular weight range' do
    # http://ctr.wikia.com/wiki/Report_how_many_SD_file_records_are_within_a_certain_molecular_weight_range
    output = wiki_code_capture_stdout do
      require 'rubabel'
      puts Rubabel.foreach("benzodiazepine.sdf.gz").count {|mol| (300..400) === mol.mol_wt }
    end
    output.should == "7\n"
    # checked with large file and it came out to 3916, which is correct
  end

  specify 'Convert SMILES file to SD file' do
    # http://ctr.wikia.com/wiki/Convert_SMILES_file_to_SD_file
    wiki_code do
      require 'rubabel'
      File.open("benzodiazepine.smi.sdf", 'w') do |out|
        Rubabel.foreach("benzodiazepine.smi.gz") do |mol|
          mol.make_3d!
          out.print mol.write_string(:sdf)
        end
      end
      # note: OpenBabel gets most stereochemistry correct.  When it cannot
      # the particular atoms are noted in a message to stderr
      # note: A minimal amount of optimization is performed in this routine
      # (following pybel's make3D routine) so that the molecule doesn't make
      # chemists reel.
    end
    unlink "benzodiazepine.smi.sdf"
  end

  specify 'Report the similarity between two structures' do
    # http://ctr.wikia.com/wiki/Report_the_similarity_between_two_structures
    output = wiki_code_capture_stdout do
      require 'rubabel'
      (mol1, mol2) = %w{CC(C)C=CCCCCC(=O)NCc1ccc(c(c1)OC)O COC1=C(C=CC(=C1)C=O)O}.map{|smile| Rubabel[smile]}
      puts mol1.tanimoto(mol2)
    end
    output.should == "0.36046511627906974\n"
    # notes: reports similarity of 0.36046511627906974
  end

  specify 'Find the 10 nearest neighbors in a data set' do
    # http://ctr.wikia.com/wiki/Find_the_10_nearest_neighbors_in_a_data_set

    output = wiki_code_capture_stdout do
      require 'rubabel'

      prev = nil
      comp = []
      Rubabel.foreach("benzodiazepine.sdf.gz").map do |mol|
        prev = mol if prev.nil?
        comp << [prev.tanimoto(mol), mol.title]
      end
      puts comp.sort.reverse[1,11].map {|v| "#{v[0].round(3)} #{v[1]}" }
    end
    output.should == "0.979 3016\n0.979 2997\n0.909 3369\n0.874 2809\n0.769 3299\n0.74 2802\n0.379 3261\n0.377 2118\n0.368 1963\n"
    # output when run on entire benzodiazepine.sdf.gz:
    #1.0 450820
    #1.0 1688
    #0.993 20351792
    #0.986 9862446
    #0.979 398658
    #0.979 398657
    #0.979 6452650
    #0.979 450830
    #0.979 44353759
    #0.979 3016
    #0.979 2997
  end

# jtp original implementation
=begin
    output = wiki_code_capture_stdout do
      require 'rubabel'
      fp1 = nil
      comparisons = Rubabel.foreach("benzodiazepine.sdf.gz").map do |mol|
        if fp1
          [OpenBabel::OBFingerprint.tanimoto(fp1, mol.ob_fingerprint), mol.title]
        else
          (fp1 = mol.ob_fingerprint) && nil
        end
      end.compact
      puts comparisons.sort.reverse[0,10].map {|v| "#{v[0].round(3)} #{v[1]}" }
    end
    output.should == "0.979 3016\n0.979 2997\n0.909 3369\n0.874 2809\n0.769 3299\n0.74 2802\n0.379 3261\n0.377 2118\n0.368 1963\n"
    # notes: still need a first class fingerprint object so this cleaner
    
    # output when run on entire benzodiazepine.sdf.gz:
    # 1.0 623918
    # 1.0 450820
    # 0.993 20351792
    # 0.986 9862446
    # 0.979 398658
    # 0.979 398657
    # 0.979 6452650
    # 0.979 450830
    # 0.979 44353759
    # 0.979 3016

    # this iterator method gives a segfault!! (not sure why)
    #mol_iter = Rubabel.foreach("benzodiazepine.sdf.gz")
    #fp1 = mol_iter.next.ob_fingerprint
    #loop do
    #  mol=mol_iter.next rescue break
    #  comparisons << [OpenBabel::OBFingerprint.tanimoto(fp1, mol.ob_fingerprint).round(2), mol.title]
    #end
    #puts comparisons.sort.reverse.map {|pair| pair.join(" ") }
  end
=end

  specify 'Depict a compound as an image' do
    # http://ctr.wikia.com/wiki/Depict_a_compound_as_an_image
    wiki_code do
      require 'rubabel'
      mol = Rubabel["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
      mol.draw(:title=>"Caffeine", :size => 900)
      mol.draw(:title=>"Caffeine", :size => 300, format: 'png', filename: 'mol.png')

      #Rubabel["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"].write_file("caffeine.svg")
    end
    #OR, using commandline
		#%x{obabel -:"CN1C=NC2=C1C(=O)N(C(=O)N2C)C" -O "mol.png" -xP 300}
  end

  xspecify 'Highlight a substructure in the depiction' do
    smarts = "c1ccc2c(c1)C(=NCCN2)c3ccccc3"
    mol = Rubabel.foreach("benzodiazepine.sdf.gz").find {|mol| mol.title == "3016" }
    p OpenBabel::OBConversion.constants
    p mol.obconv.get_options(OpenBabel::OBConversion.(:ALL))
    #mol.write("hello.svg")
  end


end


=begin
  p OpenBabel::OBOp.methods - Object.methods
  pgen = OpenBabel::OBOp.find_type("gen3D")
  p pgen.methods - Object.new.methods
  pgen.do(mol.ob)
  out.print mol.write_string(:sdf)
=end
