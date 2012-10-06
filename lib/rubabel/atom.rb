require 'matrix'
require 'andand'

require 'rubabel/bond'

class OpenBabel::OBAtom
  def upcast
    Rubabel::Atom.new(self)
  end
end

module Rubabel
  class Atom
    include Enumerable

    class << self
      # takes an element symbol and creates that atom.  If el_sym is set to
      # nil or 0, then an atom of atomic number 0 is used
      def [](el_sym=:h, id=nil)
        ob_atom = OpenBabel::OBAtom.new
        ob_atom.set_id(id) if id
        ob_atom.set_atomic_num(Rubabel::EL_TO_NUM[el_sym] || 0)
        self.new(ob_atom)
      end
    end

    # returns the molecule that is parent of this atom
    def mol
      @ob.get_parent.andand.upcast
    end

    # the OpenBabel::OBAtom object
    attr_accessor :ob

    def initialize(obatom)
      @ob = obatom
    end

    def id
      @ob.get_id
    end

    def id=(val)
      @ob.set_id(val)
    end

    # index of the atom (begins with 1)
    def idx
      @ob.get_idx
    end

    # abbreviated name, all lowercase as a Symbol
    def el
      NUM_TO_EL[atomic_num]
    end

    # abbreviated name, properly capitalized and as a String
    def element
      NUM_TO_ELEMENT[atomic_num]
    end

    # creates a bond and adds it to both atoms
    def add_atom!(other)
      obbond = OpenBabel::OBBond.new
      obbond.set_begin(self.ob)
      obbond.set_end(other.ob)
      @ob.add_bond(obbond)
      other.ob.add_bond(obbond)
      self
    end

    def each_bond(&block)
      block or return enum_for(__method__)
      iter = @ob.begin_bonds
      _bond = @ob.begin_bond(iter)
      while _bond
        block.call _bond.upcast
        _bond = @ob.next_bond(iter)
      end
    end


    # retrieves the bond 
    def get_bond(atom)
      @ob.get_bond(atom.ob).andand.upcast
    end

    # returns the bonds.  Consider using each_bond.
    def bonds
      each_bond.map.to_a
    end

    # iterates through each neighboring atom
    def each_atom(&block)
      block or return enum_for(__method__)
      iter = @ob.begin_bonds
      _atom = @ob.begin_nbr_atom(iter)
      while _atom
        block.call _atom.upcast
        _atom = @ob.next_nbr_atom(iter)
      end
    end
    alias_method :each, :each_atom

    # returns the neighboring atoms.  Consider using each_atom.
    def atoms
      each_atom.map.to_a
    end

    def atomic_mass
      @ob.get_atomic_mass
    end

    def atomic_num
      @ob.get_atomic_num
    end

    def exact_mass
      @ob.get_exact_mass
    end

    def formal_charge
      @ob.get_formal_charge
    end
    alias_method :charge, :formal_charge

    def formal_charge=(val)
      @ob.set_formal_charge(val)
    end
    alias_method :charge=, :formal_charge=

    def heavy_valence
      @ob.get_heavy_valence
    end

    def hetero_valence
      @ob.get_hetero_valence
    end

    def hyb
      @ob.get_hybridization
    end

    def implicit_valence
      @ob.get_implicit_valence
    end

    def isotope
      @ob.get_isotope
    end

    def partial_charge
      @ob.get_partial_charge
    end

    # permanently removes a hydrogen by properly incrementing the
    # spin_multiplicity (and deleting a hydrogen if one is explicitly attached
    # to the atom).  If called twice a carbene or nitrene can be made (giving
    # a spin_multiplicity of 3)
    def remove_an_h!
      new_spin = 
        case @ob.get_spin_multiplicity
        when 0 then 2
        when 2 then 3
        end
      @ob.set_spin_multiplicity(new_spin)
      mol.title = mol.to_s
      self.mol.write("before_REM.svg")
      # some molecules do not have explicit hydrogens, but after changing spin
      # they do have explicit hydrogens on *that* atom!
      atoms.each do |atom|
        if atom.atomic_num == 1
          self.mol.delete_atom(atom)
          break
        end
      end
      mol.title = mol.to_s
      self.mol.write("after_REM.svg")
      self
    end

    def ==(other)
      mol == other.mol && id == other.id 
    end

    # opposite of remove_an_h!
    # THIS IS STILL BROKEN!!!
    def add_an_h!
      abort 'bad method right now'
      new_spin = 
        case @ob.get_spin_multiplicity
        when 2 then 0
        when [1,3] then 2
        end
      h = mol.add_atom!(1)
      h.ob.set_type 'H'
      mol.add_bond!(self, h)

      @ob.set_spin_multiplicity(new_spin)
      self
    end

    def spin
      @ob.get_spin_multiplicity
    end

    def type
      @ob.get_type
    end

    def valence
      @ob.get_valence
    end

    def vector
      @ob.get_vector
    end

    def hydrogen?() @ob.is_hydrogen end
    def carbon?() @ob.is_carbon end
    def nitrogen?() @ob.is_nitrogen end
    def oxygen?() @ob.is_oxygen end
    def sulfur?() @ob.is_sulfur end
    def phosphorus?() @ob.is_phosphorus end
    def aromatic?() @ob.is_aromatic end
    def in_ring?() @ob.is_in_ring end
    def in_ring_size?() @ob.is_in_ring_size end
    def heteroatom?() @ob.is_heteroatom end
    def not_c_or_h?() @ob.is_not_cor_h end
    def connected?() @ob.is_connected end
    def one_three?() @ob.is_one_three end
    def one_four?() @ob.is_one_four end
    def carboxyl_oxygen?() @ob.is_carboxyl_oxygen end
    def phosphate_oxygen?() @ob.is_phosphate_oxygen end
    def sulfate_oxygen?() @ob.is_sulfate_oxygen end
    def nitro_oxygen?() @ob.is_nitro_oxygen end
    def amide_nitrogen?() @ob.is_amide_nitrogen end
    def polar_hydrogen?() @ob.is_polar_hydrogen end
    def non_polar_hydrogen?() @ob.is_non_polar_hydrogen end
    def aromatic_noxide?() @ob.is_aromatic_noxide end
    def chiral?() @ob.is_chiral end
    def axial?() @ob.is_axial end
    def clockwise?() @ob.is_clockwise end
    def anti_clockwise?() @ob.is_anti_clockwise end
    def positive_stereo?() @ob.is_positive_stereo end
    def negative_stereo?() @ob.is_negative_stereo end
    def chirality_specified?() @ob.has_chirality_specified end
    def chiral_volume?() @ob.has_chiral_volume end
    def hbond_acceptor?() @ob.is_hbond_acceptor end
    def hbond_donor?() @ob.is_hbond_donor end
    def hbond_donor_h?() @ob.is_hbond_donor_h end

    def double_bond?
      each_bond.any? {|bond| bond.bond_order == 2 }
    end

    def single_bond?
      each_bond.any? {|bond| bond.bond_order == 1 }
    end

    def carboxyl_carbon?
      each_atom.any?(&:carboxyl_oxygen?)
    end

    def carbonyl_oxygen?
      ats = atoms
      ats.size == 1 && ats.first.el == :c && double_bond?
    end

    def carbonyl_carbon?
      each_atom.any?(&:carbonyl_oxygen?)
    end

#    # does this carbon hold a primary alcohol
#    def primary_alcohol_carbon?
#    end

    def coords
      Vector[@ob.x, @ob.y, @ob.z]
    end

    def inspect
      "<#{type} id:#{id}>"
    end
  end
end
