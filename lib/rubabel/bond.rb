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

    attr_accessor :obbond

    def initialize(obbond)
      @obbond = obbond
    end

    # does it include the Rubabel::Atom (requires that it be holding the
    # identical obatom, i.e. uses #equal?)
    def include?(atom)
      ob_include?(atom.obatom)
    end

    def each_atom(&block)
      block or return enum_for(__method__)
      block.call @obbond.get_begin_atom.upcast
      block.call @obbond.get_end_atom.upcast
      self
    end

    alias_method :each, :each_atom

    # returns an array of Rubabel::Atoms
    def atoms
      [@obbond.get_begin_atom.upcast, @obbond.get_end_atom.upcast]
    end

    # does it include the OpenBabel::OBAtom (requires identity to be equal,
    # i.e., #equal?)
    def ob_include?(obatom)
      @obbond.get_begin_atom.equal?(obatom) || @obbond.get_end_atom.equal?(obatom)
    end

    def inspect
      "[#{atoms.map(&:inspect).join('-')}]"
    end

    # returns an array of molecules after deleting the given bond
    def split(*bonds)
    end

  end
end
