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

    attr_accessor :obatom

    def initialize(obatom)
      @obatom = obatom
    end

    def id
      @obatom.get_id
    end

    def id=(val)
      @obatom.set_id(val)
    end

    def each_bond(&block)
      block or return enum_for(__method__)
      iter = @obatom.begin_bonds
      while (_bond = @obatom.next_bond(iter))
        block.call _bond.upcast
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
      iter = @obatom.begin_bonds
      _atom = @obatom.begin_nbr_atom(iter)
      while _atom
        block.call _atom.upcast
        _atom = @obatom.next_nbr_atom(iter)
      end
    end

    # returns the neighboring atoms.  Consider using each_atom.
    def atoms
      each_atom.map.to_a
    end

    def atomic_mass
      @obatom.get_atomic_mass
    end

    def atomic_num
      @obatom.get_atomic_num
    end

    def exact_mass
      @obatom.get_exact_mass
    end

    def formal_charge
      @obatom.get_formal_charge
    end
    alias_method :charge, :formal_charge

    def heavy_valence
      @obatom.get_heavy_valence
    end

    def hetero_valence
      @obatom.get_hetero_valence
    end

    def hyb
      @obatom.get_hybridization
    end


    # returns the next Rubabel::Atom in the molecule
    def next_atom
      @obatom.get_next_atom.upcast
    end

    # index of the atom (begins with 1)
    def idx
      @obatom.get_idx
    end

    def id
      @obatom.get_id
    end

    def implicit_valence
      @obatom.get_implicit_valence
    end

    def isotope
      @obatom.get_isotope
    end

    def partial_charge
      @obatom.get_partial_charge
    end

    def spin
      @obatom.get_spin_multiplicity
    end

    def type
      @obatom.get_type
    end

    def valence
      @obatom.get_valence
    end

    def vector
      @obatom.get_vector
    end

    #def coords
    #end

    def coords
      # would like to implement with get_coordinate
      Vector[@obatom.get_x, @obatom.get_y, @obatom.get_z]
    end

    def method_missing(methd, *args, &block)
      if @obatom.respond_to?(methd)
        @obatom.send(methd, *args, &block)
      else
        super(methd, *args, &block)
      end
    end


  end
end
