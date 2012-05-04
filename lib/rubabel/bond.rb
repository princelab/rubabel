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

    # returns an array of Rubabel::Atoms
    def atoms
      [@ob.get_begin_atom.upcast, @ob.get_end_atom.upcast]
    end

    def inspect
      "[#{atoms.map(&:inspect).join('-')}]"
    end

  end
end
