require 'rubabel/atom'

class OpenBabel::OBBond
  def upcast
    Rubabel::Bond.new(self)
  end
end


module Rubabel

  # delegates to obbond object if the method is missing
  class Bond
    include Enumerable

    class << self
=begin
      def [](atom1, atom2, bond_order=1, index=0)
        abort 'cannot get bond generated properly yet, arghhh!'
        flags = 0
        obbond = OpenBabel::OBBond.new
        p obbond
        obbond.set(index, atom1.ob, atom2.ob, bond_order, flags)
        #obbond.set_end(atom2.ob)
        p obbond
        puts "GEGIN:"
        p obbond.get_begin_atom.get_id
        puts "END:"
        p obbond.get_end_atom.get_id
        obbond.set_bond_order(bond_order)
        self.new(obbond)
      end
=end
    end

    attr_accessor :ob

    def initialize(obbond)
      @ob = obbond
    end

    # considered included if the atom ids match
    def include?(atom)
      # atoms.any? {|atm| atom.id == atm.id }
      (@ob.get_begin_atom.get_id == atom.id) || (@ob.get_end_atom.get_id == atom.id)
    end

    def each_atom(&block)
      block or return enum_for(__method__)
      block.call @ob.get_begin_atom.upcast
      block.call @ob.get_end_atom.upcast
      self
    end

    alias_method :each, :each_atom

    def bond_order
      @ob.get_bond_order
    end

    # 1 = single, 2 = double, 5 = aromatic
    def bond_order=(val=1)
      @ob.set_bond_order(val)
    end

    # returns an array of Rubabel::Atoms
    def atoms
      [@ob.get_begin_atom.upcast, @ob.get_end_atom.upcast]
    end

    def inspect
      "[#{atoms.map(&:inspect).join('-')}]"
    end

  end
end
