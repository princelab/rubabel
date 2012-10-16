#encoding: utf-8

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

    # returns self
    def set_atoms!(beg_atom, end_atom)
      @ob.set_begin(beg_atom.ob)
      @ob.set_end(end_atom.ob)
      self
    end

    # Sets the beginning atom of the bond to atom. returns self
    def set_begin!(atom)
      @ob.set_begin(atom.ob)
      self
    end

    # Sets the end atom of the bond to the given atom. returns self
    def set_end!(atom)
      @ob.set_end(atom.ob)
      self
    end

    # returns an array of Rubabel::Atoms
    def atoms
      [@ob.get_begin_atom.upcast, @ob.get_end_atom.upcast]
    end

    def inspect
      bond_symbol = case bond_order
      when 2 then '='
      when 3 then '≡'
      else 
        '-'
      end
      "#{atoms.map(&:inspect).join(bond_symbol)}"
    end

    # returns self
    def +(val)
      # do we need to check the bounds here?
      newval = @ob.get_bond_order + val
      @ob.set_bond_order(newval)
      self
    end

    # won't decrease below zero. returns self
    def -(val)
      newval = @ob.get_bond_order - val
      newval = 0 if newval < 0
      @ob.set_bond_order(newval)
      self
    end

  end
end
