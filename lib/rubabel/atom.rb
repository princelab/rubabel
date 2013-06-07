#encoding: utf-8

require 'matrix'
require 'andand'

require 'rubabel/bond'

class OpenBabel::OBAtom
  def upcast
    Rubabel::Atom.new(self)
  end
end

module Rubabel
  class Atom
    include Enumerable

    class << self
      # I cannot figure out how to get the atom attached properly into a
      # molecule.  Until I figure that out, no sense in having this method
      # exposed:
      #
      # takes an atomic number or element symbol and creates that atom.  If
      # arg is set to nil, then an atom of atomic number 0 is used
      # (represented by '*' in smiles)
      #def [](arg=:c)
      #  ob_atom = OpenBabel::OBAtom.new
      #  atomic_number = 
      #    if arg.nil? then 0
      #    elsif arg.is_a?(Symbol) 
      #      Rubabel::EL_TO_NUM[arg]
      #    else
      #      arg
      #    end
      #  ob_atom.set_atomic_num(atomic_number)
      #  self.new(ob_atom)
      #end
    end

    # returns the molecule that is parent of this atom
    def mol
      @ob.get_parent.andand.upcast
    end

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

    # elemental symbol, properly capitalized and returned as a Symbol
    def element
      NUM_TO_ELEMENT[atomic_num]
    end
    alias_method :el, :element

    # connects the atom-like specifier to this atom through Molecule#add_atom!
    # returns the atom that was just added for chaining.  Takes any argument
    # that Molecule#add_atom! will take.
    def bond!(arg, bond_order=1)
      mol.add_atom!(arg, bond_order, self)
    end
    alias_method :<<, :bond!

    # connects a Rubabel::Atom object with a bond
    def connect!(atom, bond_order=1)
      @ob.get_parent.add_bond(@ob.get_idx, atom.ob.get_idx, bond_order)
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


    # retrieves the bond 
    def get_bond(atom)
      @ob.get_bond(atom.ob).andand.upcast
    end

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
    alias_method :each, :each_atom

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

    # not recognizing get_heavy_valence right now for some reason
    #def heavy_valence
    #  @ob.get_heavy_valence
    #end

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

    # doesn't take into account pH (use Molecule#do_without_hydrogens)
    def do_without_hydrogens(&block)
      _obmol = @ob.get_parent
      had_hydrogens = _obmol.has_hydrogens_added
      _obmol.delete_hydrogens(self.ob) if had_hydrogens
      reply = block.call(had_hydrogens)
      _obmol.add_hydrogens(self.ob) if had_hydrogens
      reply
    end

    # doesn't take into account pH (use Molecule#do_with_hydrogens)
    def do_with_hydrogens(&block)
      _obmol = @ob.get_parent
      had_hydrogens = _obmol.has_hydrogens_added
      _obmol.add_hydrogens(self.ob) unless had_hydrogens
      reply = block.call(had_hydrogens)
      _obmol.delete_hydrogens(self.ob) unless had_hydrogens
      reply
    end

    def remove_a_proton!(add_placeholder_hydrogens=false)
      remove_a_hydrogen!(0, add_placeholder_hydrogens)
    end

    def remove_a_hydrogen!(with_num_electrons=1, add_placeholder_hydrogens=false)
      self.dec_implicit_valence!
      case with_num_electrons
      when 0
        self.charge -= 1
      when 1
        raise NotImplementedError, "not doing free radicals just yet"
      when 2
        self.charge += 1
      end
      if @ob.explicit_hydrogen_count > 0
        _obmol = @ob.get_parent
        each do |atom|
          if atom.hydrogen?
            _obmol.delete_atom(atom.ob, false)
            break
          end
        end
      else
        if add_placeholder_hydrogens
          @ob.get_parent.add_hydrogens(@ob)
        end
      end
      self
    end

    # removes a proton with its electrons from an atom.
    # This gives precisely the same molecule as if the molecule were input by
    # smiles.  csmiles, formula, exact_mass, valence, implicit_valence, etc.
    #
    #     mol = Rubabel["CC"]
    #     mol[1].remove_a_hydride!
    #     mol == Rubabel["C[CH2+]"]  # in all characteristics
    #
    # Note, however, that with explicit hydrogens, the <b>partial charge</b>
    # is not matching up, even though every other property seems to be.  I'm
    # not sure why this is.
    #
    #     mol == Rubabel["C[CH2+]"].add_h!
    #     mol[1].remove_a_hydride!
    #     mol == Rubabel["C[CH2+]"].add_h!  # identical except in partial charge!
    def remove_a_hydride!(add_placeholder_hydrogens=false)
      remove_a_hydrogen!(2, add_placeholder_hydrogens)
    end

    # philosophy on equality: there are *so* many ways for two atoms to be
    # different that we can never really ensure that "equivalence" is met
    # without calling ~20 methods.  We narrowly define equivalence so it is
    # useful for that case and let the user make more complicated
    # equivalency/equality definitions themselves.

    # the exact same atom in the same molecule.  The equivalency test for
    # molecules is a little pricey, so better to use something like atom.id ==
    # other.id if you know you are working within the same molecule.
    def equal?(other)
      other.respond_to?(:mol) && mol.equal?(other.mol) && id == other.id 
    end

    alias_method :==, :equal?
    alias_method :eql?, :equal?

    ## opposite of remove_an_h!
    ## THIS IS STILL BROKEN!!!
    # maybe need to change the type?? C+ -> C2 or C3, but this gets really
    # invasive...  Why is this so flippin hard to do!!
    #def add_an_h!(remove_charge=true)
      #new_spin = 
        #case @ob.get_spin_multiplicity
        #when 2 then 0
        #when [1,3] then 2
        #end
      #@ob.set_spin_multiplicity(new_spin)
      
      #puts self.inspect_internals
      #puts "EXAMIN B:"
      #p self
      #p self.charge
      #(self.charge = self.charge - 1) if remove_charge
      #puts "EXAMIN A:"
      #puts self.inspect_internals
      #p self
      #p self.charge
      #puts "BEFORE:"
      #p mol.formula
      #p mol.atoms
      #mol.add_bond!(self, mol.add_atom!(1))
      #puts "AFTER:"
      #p mol.formula
      #p mol.atoms
      #abort 'here'
      #self
    #end

    def spin
      @ob.get_spin_multiplicity
    end

    def spin=(val)
      @ob.set_spin_multiplicity(val)
    end

    def type
      @ob.get_type
    end

    # Returns the current number of explicit connections.  Don't confuse this for
    # implicit_valence.
    def valence
      @ob.get_valence
    end

    # maximum number of connections expected for this atom
    def implicit_valence
      @ob.get_implicit_valence
    end

    # set the maximum number of connections expected for this atom
    def implicit_valence=(val)
      @ob.set_implicit_valence(val)
    end

    # increase by one the maximum number of connections expected for this atom
    def inc_implicit_valence!
      @ob.increment_implicit_valence
    end

    # decrease by one the maximum number of connections expected for this atom
    def dec_implicit_valence!
      @ob.decrement_implicit_valence
    end

    def vector
      @ob.get_vector
    end

    def hydrogen?() @ob.is_hydrogen end
    def carbon?() @ob.is_carbon end
    def nitrogen?() @ob.is_nitrogen end
    def oxygen?() @ob.is_oxygen end
    def sulfur?() @ob.is_sulfur end
    def phosphorus?() @ob.is_phosphorus end
    def aromatic?() @ob.is_aromatic end
    def in_ring?() @ob.is_in_ring end
    def in_ring_size?() @ob.is_in_ring_size end
    def heteroatom?() @ob.is_heteroatom end
    def not_c_or_h?() @ob.is_not_cor_h end
    def connected?() @ob.is_connected end
    def one_three?() @ob.is_one_three end
    def one_four?() @ob.is_one_four end
    def carboxyl_oxygen?() @ob.is_carboxyl_oxygen end
    def phosphate_oxygen?() @ob.is_phosphate_oxygen end
    def sulfate_oxygen?() @ob.is_sulfate_oxygen end
    def nitro_oxygen?() @ob.is_nitro_oxygen end
    def amide_nitrogen?() @ob.is_amide_nitrogen end
    def polar_hydrogen?() @ob.is_polar_hydrogen end
    def non_polar_hydrogen?() @ob.is_non_polar_hydrogen end
    def aromatic_noxide?() @ob.is_aromatic_noxide end
    def chiral?() @ob.is_chiral end
    def axial?() @ob.is_axial end
    def clockwise?() @ob.is_clockwise end
    def anti_clockwise?() @ob.is_anti_clockwise end
    def positive_stereo?() @ob.is_positive_stereo end
    def negative_stereo?() @ob.is_negative_stereo end
    def chirality_specified?() @ob.has_chirality_specified end
    def chiral_volume?() @ob.has_chiral_volume end
    def hbond_acceptor?() @ob.is_hbond_acceptor end
    def hbond_donor?() @ob.is_hbond_donor end
    def hbond_donor_h?() @ob.is_hbond_donor_h end

    # the total number of hydrogens bonded to the atom (implicit + explicit)
    def hydrogen_count
      @ob.implicit_hydrogen_count + @ob.explicit_hydrogen_count
    end

    def implicit_hydrogen_count
      @ob.implicit_hydrogen_count
    end

    def explicit_hydrogen_count
      @ob.explicit_hydrogen_count
    end
    alias_method :num_h, :hydrogen_count

    def double_bond?
      each_bond.any? {|bond| bond.bond_order == 2 }
    end

    def single_bond?
      each_bond.any? {|bond| bond.bond_order == 1 }
    end

    def carboxyl_carbon?
      each_atom.any?(&:carboxyl_oxygen?)
    end

    def carbonyl_oxygen?
      ats = atoms
      ats.size == 1 && ats.first.el == :C && double_bond?
    end

    def carbonyl_carbon?
      each_atom.any?(&:carbonyl_oxygen?)
    end

#    # does this carbon hold a primary alcohol
#    def primary_alcohol_carbon?
#    end

    def coords
      Vector[@ob.x, @ob.y, @ob.z]
    end

    def inspect
      "<#{type} id:#{id}>"
    end

    def inspect_internals
      "<" << @ob.methods.grep(/get_/).map do |mthd| 
        begin
          "#{mthd.to_s.sub(/get_/,'')}=#{@ob.send(mthd)}" 
        rescue ArgumentError
          nil
        end
      end.compact.join(" ") << ">"
    end

  end
end
