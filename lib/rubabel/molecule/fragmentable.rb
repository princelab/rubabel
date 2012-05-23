
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

      # to ensure proper fragmentation, will add_h!(ph) first at the given ph
      # an empty array is returned if there are no fragments generated.
      def fragment(opts={})
        opts = DEFAULT_OPTIONS.merge(opts)

        mol = self.dup
        mol.correct_for_ph!(opts[:ph])
        mol.remove_h!

        p 'wrd'
        rules = opts[:rules]
        fragments = []
        mol.each_match("C(O)").each do |_atoms|
        p 'hel'
          carbon = _atoms[0]
          oxygen = _atoms[1]
          # water loss
          if rules.include?(:h2oloss)
            # should there be possible water loss off of carboxylic acid???
            p mol
            p(  carbon.each_atom.map {|at| at.ob.is_carboxyl_oxygen }.count(true).size >= 2 )
            p carbon.each_atom.map {|at| at.ob.is_carboxyl_oxygen }.count(true)

            unless carbon.each_atom.map {|at| at.ob.is_carboxyl_oxygen }.count(true).size >= 2
              puts "HAFASDF"
              p carbon.get_bond(oxygen)
              frags = mol.split(carbon.get_bond(oxygen))
              p frags
              abort 'herewsdfsdfsdf'
              fragments.push *frags
            end
          end
        end
        p fragments
        fragments

        fragments.each(&:add_h!) if self.h_added?
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

