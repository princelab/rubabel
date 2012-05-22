
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

        self.correct_for_ph!(opts[:ph])

        has_hydrogens_added = h_added?
        self.remove_h! if has_hydrogens_added


        rules = opts[:rules]
        fragments = []
        if rules.include?(:co)
          self.each_match("[C,N,P]C(O)[C,N,P]").flat_map do |_atoms|
            carbon = _atoms[1]
            oxygen = _atoms[2]

            # the oxy_bonds should only have the one oxygen in it
            oxy_bonds, non_oxy_bonds = carbon.each_bond.partition {|bond| bond.include?(oxygen) }

            non_oxy_bonds.each do |bond|
              @ob.delete_bond(bond.ob)
            end





            frags = non_oxy_bonds.flat_map do |bond| 
              self.split(bond) do |fragments|
                fragments
              end
            end
             fragments.push *frags

            #oxy_bonds.first.bond_order = 2


            #fragments.push *non_oxy_bonds.flat_map {|bond| self.split(bond) }
          end
        end

        if has_hydrogens_added
          fragments.each(&:add_h!) 
          self.add_h!
        end
        fragments
      end
    end

    include Fragmentable
  end
end


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
