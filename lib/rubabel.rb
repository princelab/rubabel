require 'rubabel/version'
require 'openbabel'
require 'rubabel/molecule'

%w(atom molecule fingerprint smarts molecule_data).each do |klass|
  require "rubabel/#{klass}"
end

module Rubabel
  # the mass of an electron
  MASS_E = 0.0005486
  # www.mikeblaber.org/oldwine/chm1045/notes/Atoms/.../Atoms03.htm

  # available force-fields (would like to generate this with introspection)
  AVAILABLE_FORCEFIELDS = [:mmff94, :ghemical, :mm2, :uff] 
  DEFAULT_FORCEFIELD = AVAILABLE_FORCEFIELDS.first

  BUILDER = OpenBabel::OBBuilder.new

  # the command to execute the utility.  They are initialized to be eponymous.
  CMD = {
    babel: 'babel',
    obabel: 'obabel',
  }

  class << self

    # accepts a string specifying the molecule (calling Rubabel::Molecule.from_string) 
    # or an id (calls Rubabel::Molecule.from_id)
    def [](string, type=Rubabel::Molecule::DEFAULT_IN_TYPE)
      methd = 
        if type && Rubabel::Molecule::ID_TYPE_KEYS.include?(type)
          :from_id
        else
          :from_string
        end
      Rubabel::Molecule.send(methd, string, type)
    end

    def force_field(type=DEFAULT_FORCEFIELD)
      OpenBabel::OBForceField.find_force_field(type.to_s)
    end

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

    # returns the formats retrievable by url lookup of the id or key
    def id_formats
      Rubabel::Molecule::ID_TYPES
    end

    # returns the format Symbol that can be used for conversion, or nil if
    # the extension is not recognized.
    def format_from_ext(filename)
      obformat = OpenBabel::OBConversion.format_from_ext(filename)
      obformat.get_id.to_sym if obformat
    end

    alias_method :filetype, :format_from_ext

    # returns a format Symbol that can be used for conversion, or nil if the
    # mime-type is not recognized
    def format_from_mime(mime_type)
      obformat = OpenBabel::OBConversion.format_from_mime(mime_type)
      obformat.get_id.to_sym if obformat
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
    def molecule_from_file(filename, type=nil)
      (obmol, obconv, not_at_end) = read_first_obmol(filename, type).first
      Rubabel::Molecule.new(obmol)
    end

    # reads one molecule from the string
    def molecule_from_string(string, type=Rubabel::Molecule::DEFAULT_IN_TYPE)
      Rubabel::Molecule.from_string(string, type)
    end

    # reads the first entry and returns the OBMol object, the OBConversion object,
    # and the boolean not_at_end.  This method is not intended for public usage
    # but is necessary based on discrepancies between accessing the first
    # molecule and subsequent molecules.
    def read_first_obmol(filename, type=nil)
      type ||= format_from_ext(filename)
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
  # capitalized Symbols
  ELEMENTS = %w(H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Uut Fl Uup Lv Uus Uuo).map(&:to_sym)

  # atomic number to properly capitalized element abbreviation (as Symbol)
  NUM_TO_ELEMENT = Hash[ ELEMENTS.each_with_index.map {|el,i| [i+1,el] } ]
  # atomic number to properly capitalized element abbreviation (as Symbol)

  # the SMILES aromatic elements, listed in proper capitalized notation (e.g., :Se)
  AROMATIC_ELEMENTS = [:C, :O, :S, :Se, :N]
  # (http://esc.syrres.com/esc/docsmile.htm)

  # Along with properly capitalized element symbols (e.g., :Se) ELEMENT_TO_NUM
  # will include keys to lowercase versions of the AROMATIC_ELEMENTS
  ELEMENT_TO_NUM = NUM_TO_ELEMENT.invert

  AROMATIC_ELEMENTS.each do |el|
    ELEMENT_TO_NUM[el.to_s.downcase.to_sym] = ELEMENT_TO_NUM[el]
  end
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


