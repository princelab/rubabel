require 'openbabel'

module Rubabel
  class Molecule
    attr_accessor :obmol

    def initialize(obmol)
      @obmol = obmol
    end

    def charge
      @obmol.get_total_charge
    end

    def charge=(v)
      @obmol.set_total_charge(v)
    end

    def conformers
      @obmol.get_conformers
    end

    def add_h!
      @obmol.add_hydrogens
    end
    #alias_method :add_h!, :add_hydrogens!

    def remove_h!
      @obmol.delete_hydrogens
    end

    def separate!
      @obmol.separate
    end

    def each_atom(&block)
      # could use the native iterator in the future
      block or return enum_for(__method__)
      (1..@obmol.num_atoms).each do |n|
        block.call( @obmol.get_atom(n) )
      end
    end
    alias_method :each, :each_atom

    #def each_bond(&block)

    def cacl_desc(descnames=[])
    end

    def fingerprint(type='FP2')
    end

    alias_method :calc_fp, :fingerprint


  end
end
