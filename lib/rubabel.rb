require 'openbabel'

%w(atom molecule fingerprint smarts molecule_data).each do |klass|
  require "rubabel/#{klass}"
end

module Rubabel
  # the command to execute the utility.  They are initialized to be eponymous.
  CMD = {
    babel: 'babel',
    obabel: 'obabel',
  }

  class << self

    # returns a hash keyed by type (Symbol) pointing to a description of the
    # format
    def in_formats
      @in_formats ||= formats_to_hash(OpenBabel::OBConversion.new.get_supported_input_format)
    end

    # returns a hash keyed by type (Symbol) pointing to a description of the
    # format
    def out_formats
      @out_formats ||= formats_to_hash(OpenBabel::OBConversion.new.get_supported_output_format)
    end

    # determines the extension from filename if type is nil
    def foreach(filename, type=nil, &block)
      block or return enum_for(__method__, filename, type)
      (obmol, obconv, not_at_end) = read_first_obmol(filename, type)
      # the obmol is not valid if we are already at the end!
      while not_at_end
        block.call Rubabel::Molecule.new(obmol)
        obmol = OpenBabel::OBMol.new
        not_at_end = obconv.read(obmol)
      end
    end

    # returns a Rubabel::Molecule (the first in the file if there are
    # multiples).  See ::foreach for accessing all molecules in a file
    # determines the type from the extension if type is nil.
    def read_file(filename, type=nil)
      obmol = read_first_obmol(filename, type).first
      Rubabel::Molecule.new(obmol)
    end

    # reads one molecule from the string
    def read_string(string, type=:smi)
      obmol = OpenBabel::OBMol.new
      obconv = OpenBabel::OBConversion.new
      obconv.set_in_format(type.to_s) || raise(ArgumentError, "invalid format #{type}")
      success = obconv.read_string(obmol, string)
      Rubabel::Molecule.new(obmol)
    end

    # returns a filetype symbol based on the extension
    def filetype(filename)
      # should use the openbabel method in the future
      File.extname(filename)[1..-1].to_sym
    end

    private
    # reads the first entry and returns the OBMol object, the OBConversion object,
    # and the boolean not_at_end.  This method is not intended for public usage
    # but is necessary based on discrepancies between accessing the first
    # molecule and subsequent molecules.
    def read_first_obmol(filename, type=nil)
      type ||= filetype(filename)
      obconv = OpenBabel::OBConversion.new
      obconv.set_in_format(type.to_s) || raise(ArgumentError, "invalid format #{type}")
      obmol = OpenBabel::OBMol.new
      not_at_end = obconv.read_file(obmol, filename)
      [obmol, obconv, not_at_end]
    end

    def formats_to_hash(format_strings)
      Hash[ 
        format_strings.map do |str| 
          pair = str.split(/\s+--\s+/)
          [pair[0].to_sym, pair[1]]
        end 
      ]
    end
  end

  #protected_class_method :read_first_obmol

end

module Rubabel
  # capitalized strings
  ELEMENTS = %w(H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Uut Fl Uup Lv Uus Uuo)

  # atomic number to properly capitalized element abbreviation
  NUM_TO_ELEMENT = Hash[ ELEMENTS.each_with_index.map {|el,i| [i+1,el] } ]

  # atomic number to lowercase symbol abbreviation
  NUM_TO_EL = Hash[ ELEMENTS.each_with_index.map {|el,i| [i+1,el.downcase.to_sym] } ]
end

=begin
OBConversion conv;
OBMol mol;
bool success = conv.SetInFormat("sdf");
if(success)
  {
    bool notatend = conv.ReadFile(&mol, "myfile.sdf");
    // Do something with mol
    while(notatend)
      {
        notatend = conv.Read(&mol);
        // Do something with mol
      }
  }
=end
