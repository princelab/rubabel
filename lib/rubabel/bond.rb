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
      def [](atom1, atom2)
        obbond = OpenBabel::OBBond.new
        obbond.set_begin(atom1.ob)
        obbond.set_end(atom2.ob)
        self.new(obbond)
      end
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
