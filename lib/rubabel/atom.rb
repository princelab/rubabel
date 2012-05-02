require 'matrix'

module Rubabel
  class Atom

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

    def heavy_valence
      @obatom.get_heavy_valence
    end

    def hetero_valence
      @obatom.get_hetero_valence
    end

    def hyb
      @obatom.get_hybridization
    end

    def bonds
      ################## cast..
      @obatom.get_bonds
    end

    # returns the next Rubabel::Atom in the molecule
    def next_atom
      Rubabel::Atom.new(@obatom.get_next_atom)
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

  end
end
