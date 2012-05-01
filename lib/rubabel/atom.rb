
module Rubabel
  class Atom
    attr_accessor :obatom
    def initialize(obatom)
      @obatom = obatom
    end

    %w(atomic_mass atomic_num exact_mass).each do |meth|
      define_method(meth) do 
        @obatom.send("get_#{meth}")
      end

      define_method("#{meth}=") do |val|
        @obatom.send("set_#{meth}", val)
      end
    end

    def coords
    end

    def exact_mass
    end

    def formal_charge
    end

    def heavy_valence
    end

    def hetero_valence
    end

    def hyb
    end

    # index of the atom (begins with 1)
    def idx
    end

    def implicit_valence
    end

    def isotope
    end

    def partial_charge
    end

    def spin
    end

    def type
    end

    def valence
    end

    def vector
    end

  end
end
