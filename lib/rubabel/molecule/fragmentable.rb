
module Rubabel
  class Molecule
    module Fragmentable
      RULES = [:co]
      #ADDUCTS = [:lioh, :nh4cl, :nh4oh]

      DEFAULT_OPTIONS = {
        rules: RULES,
        #adduct: nil,
        ph: 7.4,
        # return only the set of unique fragments
        uniq: false, 
      }

      # molecules and fragments should all have hydrogens added (add_h!)
      # before calling this method
      # 
      # For instance, water loss with double bond formation is not allowable
      # for NCC(O)CC => CCC=C[NH2+], presumably because of the lone pair and
      # double bond resonance.
      #     
      def allowable_fragmentation?(frags)
        self.num_atoms == frags.map(&:num_atoms).reduce(:+)
      end

      # will turn bond into a double bond, yield the changed molecule, then
      # return the bond to the original state when the block is closed
      # returns whatever the block returned
      def feint_double_bond(bond, &block)
        orig = bond.bond_order
        bond.bond_order = 2
        reply = block.call(self)
        bond.bond_order = orig
        reply
      end

      # to ensure proper fragmentation, will add_h!(ph) first at the given ph
      # an empty array is returned if there are no fragments generated.
      def fragment(opts={})
        opts = DEFAULT_OPTIONS.merge(opts)

        had_hydrogens = self.h_added?

        self.correct_for_ph!(opts[:ph])
        self.remove_h!

        rules = opts[:rules]
        fragments = []
        self.each_match("CO").each do |_atoms|
          (carbon, oxygen) = _atoms
          carbon_nbrs = carbon.atoms.reject {|atom| atom == oxygen }

          case oxygen.bonds.size
          when 1  # an alcohol
            # water loss
            double_bondable = carbon_nbrs.select {|atm| atm.type == 'C3' }
            if rules.include?(:h2oloss) && (double_bondable.size > 0) && !carbon.carboxyl_carbon?
              frag_sets = double_bondable.map do |dbl_bondable_atom|
                frags = feint_double_bond(dbl_bondable_atom.get_bond(carbon)) do |_mol|
                  # TODO: check accuracy before completely splitting for efficiency
                  frags = _mol.split(carbon.get_bond(oxygen))
                  frags.map(&:add_h!)
                end
              end

              self.add_h!
              frag_sets.select! do |_frags| 
                self.allowable_fragmentation?(_frags)
              end
              fragments.push *frag_sets.flatten(1)
            end
            # oxygen bonded to something else (per-oxide??)
            # also could be ether situation...
          when 2  
          end
        end
        unless had_hydrogens
          fragments.each(&:remove_h!)
          self.remove_h!
        end
        fragments
      end

    end
    include Fragmentable
  end
end


  #          co_bond = carbon.get_bond(oxygen)
  #            left_to_c_bond = carbon.get_bond(left)
  #            right_to_c_bond = carbon.get_bond(right)
  #
  #            co_bond.bond_order = 2
  #
  #            [left_to_c_bond, right_to_c_bond].flat_map do |other_to_c_bond|
  #              mol.ob.delete_bond(other_to_c_bond.ob, false)
  #              pieces = mol.ob.separate.map(&:upcast)
  #              mol.ob.add_bond(other_to_c_bond.ob)
  #              pieces
  #            end


=begin

        # duplicate the molecule so we can do what we like with it
        mol = self.dup

        has_hydrogens_added = h_added?
        mol.remove_h! if has_hydrogens_added

        mol.correct_for_ph!(opts[:ph])

        rules = opts[:rules]
        fragments = []
        if rules.include?(:co)
          mol.each_match("C(O)").flat_map do |_atoms|
            carbon = _atoms.first
            non_oxygen = carbon.each_bond.reject {|bond| bond.include?(_atoms.last) }
            non_oxygen.each {|bond| p mol.split(bond) }

            fragments.push *non_oxygen.flat_map {|bond| mol.split(bond) }
          end
        end
        p fragments
        abort 'here'
        fragments.each(&:add_h!) if has_hydrogens_added
        fragments
      end
    end

=end


  #            [[left_to_c_bond, left], [right_to_c_bond, right]].flat_map do |other_to_carbony_c, other|
  #              puts "INSIDE!!!"
  #              pieces = mol.split(other_to_carbony_c)
  #              c_in_pieces = nil
  #              oxy_in_pieces = nil
  #              other_in_pieces = nil
  #              pieces.each do |piece| 
  #                piece.each_atom do |atom| 
  #                  p piece
  #                  p [atom.id, other.id]
  #                  other_in_pieces = atom if atom.id == other.id
  #                  c_in_pieces = atom if atom.id == carbon.id 
  #                  oxy_in_pieces = atom if atom.id == oxygen.id 
  #                  break if c_in_pieces && oxy_in_pieces
  #                end
  #                break if c_in_pieces && oxy_in_pieces
  #              end
  #              oxygen_bond = c_in_pieces.get_bond(oxy_in_pieces)
  #              oxygen_bond.bond_order = 2
  #
  #              puts "EXAMINE:other"
  #              p other_in_pieces.ob
  #              p other_in_pieces.mol.csmiles
  #              other_mol = other_in_pieces.mol
  #              ob_atom = other_mol.ob.new_atom
  #              ob_atom.set_atomic_num 1
  #              newbond = OpenBabel::OBBond.new 
  #              ob_atom
  #

