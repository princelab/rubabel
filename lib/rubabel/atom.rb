require 'matrix'

require 'rubabel/bond'

class OpenBabel::OBAtom
  def upcast
    Rubabel::Atom.new(self)
  end
end

module Rubabel
  class Atom
    include Enumerable

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

    def each_bond(&block)
      block or return enum_for(__method__)
      iter = @ob.begin_bonds
      _bond = @ob.begin_bond(iter)
      while _bond
        block.call _bond.upcast
        _bond = @ob.next_bond(iter)
      end
    end

    alias_method :each, :each_bond

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

    def coords
      Vector[@ob.x, @ob.y, @ob.z]
    end

    def inspect
      "<#{type} id:#{id}>"
    end
  end
end
