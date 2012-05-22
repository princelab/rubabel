require 'openbabel'
require 'rubabel'
require 'rubabel/atom'
require 'rubabel/bond'
require 'rubabel/molecule/fragmentable'

class OpenBabel::OBMol
  def upcast
    Rubabel::Molecule.new(self)
  end
end

class OpenBabelUnableToSetupForceFieldError < RuntimeError
end

module Rubabel
    # yet to implement: 
  class Molecule
    DEFAULT_FINGERPRINT = "FP2"
    include Enumerable

    # the OpenBabel::OBmol object
    attr_accessor :ob

    # the OpenBabel::OBConversion object
    attr_accessor :obconv

    class << self

      def tanimoto(mol1, mol2, type=DEFAULT_FINGERPRINT)
        OpenBabel::OBFingerprint.tanimoto(mol1.ob_fingerprint(type), mol2.ob_fingerprint(type))
      end

      def from_file(file, type=nil)
        Rubabel.molecule_from_file(file, type)
      end

      def from_string(string, type=:smi)
        Rubabel.molecule_from_string(string, type)
      end
    end

    DEFAULT_OUT_TYPE = :can

    # attributes
    def title() @ob.get_title end
    def title=(val) @ob.set_title(val) end

    def charge() @ob.get_total_charge end
    def charge=(v) @ob.set_total_charge(v) end

    def spin() @ob.get_total_spin_multiplicity end

    def mol_wt() @ob.get_mol_wt end
    alias_method :avg_mass, :mol_wt

    def exact_mass() @ob.get_exact_mass end

    # returns a string representation of the molecular formula.  Not sensitive
    # to add_h!
    def formula() @ob.get_formula end


    def initialize(obmol, obconv=nil)
      @obconv = obconv
      @ob = obmol
    end

    # returns a list of atom indices matching the patterns (corresponds to
    # the OBSmartsPattern::GetUMapList() method).  Note that the original
    # GetUMapList returns atom *numbers* (i.e., the index + 1).  This method
    # returns the zero indexed indices.
    def smarts_indices(smarts_or_string)
      pattern = smarts_or_string.is_a?(Rubabel::Smarts) ? smarts_or_string : Rubabel::Smarts.new(smarts_or_string)
      pattern.ob.match(@ob)
      pattern.ob.get_umap_list.map {|atm_indices| atm_indices.map {|i| i - 1 } }
    end

    # yields atom arrays matching the pattern.  returns an enumerator if no
    # block is given
    def each_match(smarts_or_string, &block)
      block or return enum_for(__method__, smarts_or_string)
      _atoms = self.atoms
      smarts_indices(smarts_or_string).each do |ar|
        block.call(_atoms.values_at(*ar))
      end
    end

    # returns an array of matching atom sets.  Consider using each_match.
    def matches(smarts_or_string)
      each_match(smarts_or_string).map.to_a
    end

    def matches?(smarts_or_string)
      # TODO: probably a more efficient way to do this using API
      smarts_indices(smarts_or_string).size > 0
    end

    # returns an array of OpenBabel::OBRing objects.
    def ob_sssr
      @ob.get_sssr.to_a
    end

    #def conformers
      # Currently returns an object of type
      # SWIG::TYPE_p_std__vectorT_double_p_std__allocatorT_double_p_t_t
      #vec = @ob.get_conformers
    #end
    
    # are there hydrogens added yet
    def hydrogens_added?
      @ob.has_hydrogens_added
    end
    alias_method :h_added?, :hydrogens_added?

    # returns self.  Corrects for ph if ph is not nil.  NOTE: the reversal of
    # arguments from the OpenBabel api.
    def add_h!(ph=nil, polaronly=false)
      if ph.nil?
        @ob.add_hydrogens(polaronly)
      else
        @ob.add_hydrogens(polaronly, true, ph)
      end
      self
    end

    # returns self.  If ph is nil, then #neutral! is called
    def correct_for_ph!(ph=7.4)
      ph.nil? ? neutral! : @ob.correct_for_ph(ph)
      self
    end

    # simple method to coerce the molecule into a neutral charge state.
    # It does this by removing any charge from each atom and then removing the
    # hydrogens (which will then can be added back by the user and will be
    # added back with proper valence).  If the molecule had hydrogens added it
    # will return the molecule with hydrogens added
    # returns self.
    def neutral!
      had_hydrogens = h_added?
      atoms.each {|atom| atom.charge = 0 if (atom.charge != 0) }
      remove_h!
      add_h! if had_hydrogens 
      self
    end

    # adds hydrogens 
    #def add_h_at_ph!(ph=7.4)
    #  # creates a new molecule (currently writes to smiles and uses babel
    #  # commandline to get hydrogens at a given pH; this is because no pH model
    #  # in ruby bindings currently).
#
#      ## write the file with the molecule
#      #self.write_file("tmp.smi")
#      #system "#{Rubabel::CMD[:babel]} -i smi tmp.smi -p #{ph} -o can tmp.can"
#      #Molecule.from_file("tmp.can")
#    end
#    #alias_method :add_h!, :add_hydrogens!

    def remove_h!
      @ob.delete_hydrogens
    end

    # calls separate on the OBMol object
    def separate!
      @ob.separate
    end

    # returns just the smiles string :smi (not the id)
    def smiles
      to_s(:smi)
    end

    # returns just the smiles string (not the id)
    def csmiles
      to_s(:can)
    end

    # equal if their canonical SMILES strings (including ID) are equal.  This
    # may become more rigorous in the future
    def ==(other)
      self.write_string(:can) == other.write_string(:can)
    end

    # iterates over the molecule's Rubabel::Atom objects
    def each_atom(&block)
      # could use the C++ iterator in the future
      block or return enum_for(__method__)
      iter = @ob.begin_atoms
      atom = @ob.begin_atom(iter)
      while atom
        block.call atom.upcast
        atom = @ob.next_atom(iter)
      end
    end
    alias_method :each, :each_atom

    # iterates over the molecule's Rubabel::Bond objects
    def each_bond(&block)
      # could use the C++ iterator in the future
      block or return enum_for(__method__)
      iter = @ob.begin_bonds
      obbond = @ob.begin_bond(iter)
      while obbond
        block.call obbond.upcast
        obbond = @ob.next_bond(iter)
      end
      self
    end

    # returns the array of bonds.  Consider using #each_bond
    def bonds
      each_bond.map.to_a
    end

    # returns the array of atoms. Consider using #each
    def atoms
      each_atom.map.to_a
    end

    def dim
      @ob.get_dimension
    end

    #def each_bond(&block)

    # TODO: implement
    #def calc_desc(descnames=[])
    #end

    # TODO: implement
    # list of supported descriptors. (Not Yet Implemented!)
    #def descs
    #end

    def tanimoto(other, type=DEFAULT_FINGERPRINT)
      Rubabel::Molecule.tanimoto(self, other, type)
    end

    # returns a  std::vector<unsigned int> that can be passed directly into
    # the OBFingerprint.tanimoto method
    def ob_fingerprint(type=DEFAULT_FINGERPRINT)
      fprinter = OpenBabel::OBFingerprint.find_fingerprint(type) || raise(ArgumentError, "fingerprint type not found")
      fp = OpenBabel::VectorUnsignedInt.new
      fprinter.get_fingerprint(@ob, fp) || raise("failed to get fingerprint for #{mol}")
      fp
    end

    def split(*bonds)
      bonds.each do |bond|
        unless @ob.delete_bond(bond.ob, false)
          raise "#{bond.inspect} not deleted!" 
        end
      end
      frags = @ob.separate.map(&:upcast)
      bonds.each {|bond| @ob.add_bond(bond.ob) }
      frags
    end

    alias_method :separate, :split

    # emits smiles without the trailing tab, newline, or id.  Use write_string
    # to get the default OpenBabel behavior (ie., tabs and newlines).
    def to_s(type=DEFAULT_OUT_TYPE)
      string = write_string(type)
      case type
      when :smi, :smiles, :can
        # remove name with openbabel options in the future
        string.split(/\s+/).first
      else
        string
      end
    end

    # returns a Rubabel::MoleculeData hash
    def data
      Rubabel::MoleculeData.new(@ob)
    end

    # sensitive to add_h!
    def num_atoms() @ob.num_atoms  end
    def num_bonds() @ob.num_bonds  end
    def num_hvy_atoms() @ob.num_hvy_atoms  end
    def num_residues() @ob.num_residues  end
    def num_rotors() @ob.num_rotors  end

    # If filename_or_type is a symbol, then it will return a string of that type.
    # If filename_or_type is a string, will write to the filename
    # given no args, returns a DEFAULT_OUT_TYPE string
    def write(filename_or_type=:can)
      if filename_or_type.is_a?(Symbol)
        write_string(filename_or_type)
      else
        write_file(filename_or_type)
      end
    end

    # adds hydrogens if necessary.  Performs only steepest descent
    # optimization (no rotors optimized)
    # returns self
    def local_optimize!(forcefield=DEFAULT_FORCEFIELD, steps=500)
      add_h! unless hydrogens_added?
      if dim == 3
        ff = Rubabel.force_field(forcefield.to_s)
        ff.setup(@ob) || raise(OpenBabelUnableToSetupForceFieldError)
        ff.steepest_descent(steps)  # is the default termination count 1.0e-4 (used in obgen?)
        ff.update_coordinates(@ob)
      else
        make_3d!(forcefield, steps) 
      end
      self
    end

    def global_optimize!(forcefield=DEFAULT_FORCEFIELD, steps=1000)
      if dim != 3
        # don't bother optimizing yet (steps=nil)
        make_3d!(DEFAULT_FORCEFIELD, nil)
      end
    end

    # does a bit of basic local optimization unless steps is set to nil
    # returns self
    def make_3d!(forcefield=DEFAULT_FORCEFIELD, steps=50)
      BUILDER.build(@ob)
      @ob.add_hydrogens(false, true) unless hydrogens_added?
      local_optimize!(forcefield, steps) if steps
      self
    end
    alias_method :make3d!, :make_3d!

    def write_string(type=DEFAULT_OUT_TYPE)
      @obconv ||= OpenBabel::OBConversion.new
      @obconv.set_out_format(type.to_s)
      @obconv.write_string(@ob)
    end

    # writes to the file based on the extension given.  If type is given
    # explicitly, then it is used.
    def write_file(filename, type=nil)
      type ||= Rubabel.filetype(filename)
      File.write(filename, write_string(type))
    end

    def inspect
      "#<Mol #{to_s}>"
    end

  end
end

# NOTES TO SELF:
=begin
~/tools/openbabel-rvmruby1.9.3/bin/babel -L descriptors
abonds    Number of aromatic bonds
atoms    Number of atoms
bonds    Number of bonds
cansmi    Canonical SMILES
cansmiNS    Canonical SMILES without isotopes or stereo
dbonds    Number of double bonds
formula    Chemical formula
HBA1    Number of Hydrogen Bond Acceptors 1 (JoelLib)
HBA2    Number of Hydrogen Bond Acceptors 2 (JoelLib)
HBD    Number of Hydrogen Bond Donors (JoelLib)
InChI    IUPAC InChI identifier
InChIKey    InChIKey
L5    Lipinski Rule of Five
logP    octanol/water partition coefficient
MR    molar refractivity
MW    Molecular Weight filter
nF    Number of Fluorine Atoms
s    SMARTS filter
sbonds    Number of single bonds
smarts    SMARTS filter
tbonds    Number of triple bonds
title    For comparing a molecule's title
TPSA    topological polar surface area

# why can't I find any descriptors????
>> require 'openbabel'
=> true
>> OpenBabel::OBDescriptor.find_type("logP")
=> nil
>> OpenBabel::OBDescriptor.find_type("MR")
=> nil

=end

