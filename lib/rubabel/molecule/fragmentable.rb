require 'set'
require 'rubabel/core_ext/putsv'
require 'rubabel/core_ext/enumerable'

module Rubabel
  class Molecule
    module Fragmentable

      #:sp3c_oxygen_asymmetric_far_sp3, :sp3c_nitrogen_asymmetric_far_sp3,
      #RULES = Set[ :alcohol_to_aldehyde, :peroxy_to_carboxy, :co2_loss, 
      #  :sp3c_oxygen_double_bond_far_side_sp3, :sp3c_oxygen_double_bond_far_side_sp2, :sp3c_oxygen_double_bond_water_loss, :sp3c_nitrogen_double_bond,
      #]
      #ADDUCTS = [:lioh, :nh4cl, :nh4oh]
      #CO_RULES = Set[:alcohol_to_aldehyde, :peroxy_to_carboxy, :co2_loss, 
      #  :sp3c_oxygen_double_bond_water_loss, :sp3c_oxygen_double_bond_far_side_sp2, :sp3c_oxygen_double_bond_far_side_sp3, :sp3c_oxygen_asymmetric_far_sp3
      #]

      RULES = Set[:cad_o, :cad_oo, :oxe, :oxe_phde, :ohs]

      DEFAULT_OPTIONS = {
        rules: RULES,
        #adduct: nil,
        #ph: 7.4,
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

      # add_h! to self, then selects allowable fragments
      def allowable_fragment_sets!(fragment_sets)
        self.add_h!
        fragment_sets.select do |_frags| 
          putsv "ExMAIN:"
          putsv _frags.inspect
          putsv self.allowable_fragmentation?(_frags)
          self.allowable_fragmentation?(_frags)
        end
      end

      # will turn bond into a double bond, yield the changed molecule, then
      # return the bond to the original state when the block is closed
      # returns whatever the block returned
      def feint_double_bond(bond, give_e_pair=nil, get_e_pair=nil, &block)
        orig = bond.bond_order
        bond.bond_order = 2
        reply = 
          if give_e_pair || get_e_pair
            feint_e_transfer(give_e_pair, get_e_pair, &block)
          else
            block.call(self)
          end
        bond.bond_order = orig
        reply
      end

      # warning, this method adds_h! to the calling molecule
      def electrophile_snatches_electrons(carbon, electrophile)
        self.add_h!
        frags = self.split(carbon.get_bond(electrophile))
        raise NotImplementedError
        # don't check for allowable fragments because it 
        #allowable_fragment_sets!([frag_set])
      end

      def feint_e_transfer(give_e_pair=nil, get_e_pair=nil, &block)
        if give_e_pair
          gc_orig = give_e_pair.charge
          give_e_pair.charge = gc_orig + 1
        end
        if get_e_pair
          rc_orig = get_e_pair.charge
          get_e_pair.charge = rc_orig - 1
        end

        reply = block.call(self)
        
        give_e_pair.charge = gc_orig if give_e_pair
        get_e_pair.charge = rc_orig if get_e_pair

        reply
      end

      def near_side_double_bond_break(carbon, electrophile)
        frag_sets = carbon.atoms.select {|atom| atom.carbon? && !atom.aromatic? && (atom.num_h > 0) }.map do |near_c3|
          frags = feint_double_bond(carbon.get_bond(near_c3)) do |_mol|
            frags = _mol.split(electrophile.get_bond(carbon))
          end
        end
        #allowable_fragment_sets!(frag_sets)
      end

      #def alcohol_to_aldehyde(carbon, oxygen, carbon_nbrs)
      #  # alcohol becomes a ketone and one R group is released
      #  frag_sets = carbon_nbrs.select {|atom| atom.type == 'C3' }.map do |_atom|
      #    frags = feint_double_bond(carbon.get_bond(oxygen)) do |_mol|
      #      frags = _mol.split(carbon.get_bond(_atom))
      #      frags.map(&:add_h!)
      #    end
      #  end
      #  allowable_fragment_sets!(frag_sets)
      #end

      def co2_loss(carbon, oxygen, c3_nbr)
        # carboxyl rules ...
        # neutral carbon dioxide loss with anion gain on attaching group
        # (if carbon)
        frags = feint_double_bond(carbon.get_bond(oxygen), oxygen, c3_nbr) do |_mol|
          frags = _mol.split(c3_nbr.get_bond(carbon))
          frags.map(&:add_h!)
        end
        allowable_fragment_sets!([frags])
      end

      def peroxy_to_carboxy(carbon, oxygen, carbon_nbrs, oxygen_nbr)
        if oxygen_nbr.el == :o # has a neighbor oxygen
          distal_o = oxygen_nbr
          if distal_o.bonds.size == 1  # this is a peroxy
            frag_sets = carbon_nbrs.select {|atom| atom.type == 'C3' }.map do |_atom|
              self.swap!(carbon, _atom, oxygen, distal_o)
              frags = feint_double_bond(carbon.get_bond(oxygen)) do |_mol|

                # we swapped the atoms so the bond to split off is now
                # attached to the oxygen
                frags = _mol.split(oxygen.get_bond(_atom))
                frags.map(&:add_h!)
              end
              self.swap!(carbon, distal_o, oxygen, _atom)
              frags
            end
            allowable_fragment_sets!(frag_sets)
          end
        end

      end

      # splits the molecule between the carbon and carbon_nbr, adds a double
      # bond between the carbon and oxygen, and moves whatever was on the
      # oxygen (e.g., an OH or a charge) to the carbon_nbr. Returns two new
      # molecules.
      def carbonyl_oxygen_dump(carbon, oxygen, carbon_nbr)
        appendage = oxygen.atoms.find {|a| a.el != :c }
        if oxygen.charge != 0
          ocharge = oxygen.charge
        end
        nmol = self.dup
        new_oxygen = nmol.atom(oxygen.id)
        new_carbon = nmol.atom(carbon.id)
        new_carbon_nbr = nmol.atom(carbon_nbr.id)
        new_appendage = nmol.atom(appendage.id) if appendage
        nmol.delete_bond(new_carbon.get_bond(new_carbon_nbr))
        if new_appendage
          nmol.delete_bond(new_oxygen.get_bond(new_appendage)) 
          nmol.add_bond!(new_carbon_nbr, new_appendage)
        end
        if ocharge
          new_carbon_nbr.charge += ocharge
          new_oxygen.charge -= ocharge
        end
        new_carbon.get_bond(new_oxygen).bond_order = 2
        nmol.split
      end

      # breaks the bond and gives the electrons to the oxygen
      def carbon_oxygen_esteal(carbon, oxygen)
        nmol = self.dup
        ncarbon = nmol.atom(carbon.id)
        noxygen = nmol.atom(oxygen.id)
        nmol.delete_bond(ncarbon, noxygen)
        ncarbon.remove_an_h!
        #noxygen.ob.set_spin_multiplicity  1
        noxygen.spin = 1
        noxygen.charge = -1
        nmol.split
      end

      # returns the duplicated molecule and the equivalent atoms
      def dup_molecule(atoms=[])
        nmol = self.dup
        [nmol, atoms.map {|old_atom| nmol.atom(old_atom.id) }]
      end

      # returns molecules created from splitting between the electrophile and
      # the center and where the bond order is increased between the center
      # and center_nbr
      def break_with_double_bond(electrophile, center, center_nbr)
        (nmol, (nele, ncarb, ncarb_nbr)) = self.dup_molecule([electrophile, center, center_nbr])
        nmol.delete_bond(nele, ncarb)
        ncarb_nbr.get_bond(ncarb) + 1
        nmol.split
      end

      # an empty array is returned if there are no fragments generated.
      # Hydrogens are added at a pH of 7.4, unless they have already been
      # added.
      #
      #     :rules => queryable by :include? set of rules
      #     :uniq => false
      def fragment(opts={})
        only_uniqs = true
        opts = DEFAULT_OPTIONS.merge(opts)
        opts[:rules].each do |rule| 
          raise ArgumentError, "bad rule: #{rule}" unless RULES.include?(rule)
        end

        had_hydrogens = self.h_added?
        self.correct_for_ph!(7.4) unless had_hydrogens
        self.remove_h!

        fragment_sets = []

        if opts[:rules].any? {|r| [:cad_o, :cad_oo].include?(r) }
          self.each_match("C[O;h1,O]", only_uniqs) do |carbon, oxygen|
            carbon.atoms.select {|a| a.el == :c }.each do |carbon_nbr|
              fragment_sets << carbonyl_oxygen_dump(carbon, oxygen, carbon_nbr)
            end
          end
        end
        if opts[:rules].any? {|r| [:oxe].include?(r) }
          self.each_match("C-O", only_uniqs) do |carbon, oxygen|
            fragment_sets << carbon_oxygen_esteal(carbon, oxygen)
          end
        end
        # right now implemented so that a beta hydrogen has to be availabe for
        # extraction
        if opts[:rules].any? {|r| [:ohs].include?(r) }
          self.each_match("C[C,O]-O", only_uniqs) do |beta_c, center, oxygen|
            next unless beta_c.ob.implicit_hydrogen_count > 0
            fragment_sets << break_with_double_bond(oxygen, center, beta_c)
          end
        end
        if opts[:rules].any? {|r| [:oxe_phde].include?(r) }
          self.each_match("P-O-C", only_uniqs) do |phosphate, oxygen, carbon|
            frag_set = carbon_oxygen_esteal(phosphate, oxygen)
            frag_set.map! &:convert_dative_bonds!
            fragment_sets << frag_set
          end
        end
      
        unless had_hydrogens
          fragment_sets.each {|set| set.each(&:remove_h!) }
          self.remove_h!
        end
        if opts[:uniq]
          # TODO: impelent properly
          raise NotImplementedError
          #fragment_sets = fragment_sets.uniq_by(&:csmiles)
        end

        fragment_sets
      end


      #        had_hydrogens = self.h_added?

      #self.correct_for_ph!(opts[:ph])
      #self.remove_h!

      #rules = opts[:rules]
      #fragment_sets = []
      #if rules.any? {|rule| CO_RULES.include?(rule) }
      #putsv "matching C-O"
      #self.each_match("CO").each do |_atoms|
      ## note: this will *not* match C=O
      #(carbon, oxygen) = _atoms
      #carbon_nbrs = carbon.atoms.reject {|atom| atom == oxygen }
      #c3_nbrs = carbon_nbrs.select {|atm| atm.type == 'C3' }
      ## pulling this out here causes it to work incorrectly internally
      ## (not sure why)
      ##co_bond = carbon.get_bond(oxygen)

      #case oxygen.bonds.size # non-hydrogen bonds
      #when 1  # *must* be an alcohol or a carboxylic acid
      #putsv "#{csmiles} oxygen has no other bonds besides C-O (alcohol or carboxylic acid)"
      #if carbon.type == 'C3'
      #if rules.include?(:sp3c_oxygen_double_bond_water_loss) 
      #putsv "rule :sp3c_oxygen_double_bond_water_loss"
      #fragment_sets.push *near_side_double_bond_break(carbon, oxygen)
      #end
      #if rules.include?(:alcohol_to_aldehyde)
      #putsv "rule :alcohol_to_aldehyde"
      #fragment_sets.push *alcohol_to_aldehyde(carbon, oxygen, carbon_nbrs)
      #end
      #elsif carbon.carboxyl_carbon?
      #if rules.include?(:co2_loss)
      #putsv "rule :co2_loss"
      #if c3_nbr = c3_nbrs.first
      #fragment_sets.push *co2_loss(carbon, oxygen, c3_nbr)
      #end
      #end
      #end
      #when 2
      #putsv "#{csmiles} c-o & oxygen has 2 non-hydrogen bonds"
      #oxygen_nbr = oxygen.atoms.reject {|atom| atom.idx == carbon.idx }.first
      #if carbon.type == 'C3'
      #if rules.include?(:peroxy_to_carboxy)
      #fragment_sets.push *peroxy_to_carboxy(carbon, oxygen, carbon_nbrs, oxygen_nbr)
      #end
      ## ester and ethers (look *only* on close side for places to make
      ## double bond)

      #if oxygen_nbr.type == 'C3'
      #putsv "oxygen nbr is C3"
      #if rules.include?(:sp3c_oxygen_double_bond_far_side_sp3) 
      #putsv "rule :sp3c_oxygen_double_bond_far_side_sp3"
      #fragment_sets.push *near_side_double_bond_break(carbon, oxygen)
      #end
      #if rules.include?(:sp3c_oxygen_asymmetric_far_sp3)
      #putsv "rule :sp3c_oxygen_asymmetric_far_sp3"
      ## only returns a single frag set
      #fragment_sets.push electrophile_snatches_electrons(carbon, oxygen)
      #end
      #end
      #if oxygen_nbr.type == 'C2'
      #if rules.include?(:sp3c_oxygen_double_bond_far_side_sp2)
      #putsv "rule :sp3c_oxygen_double_bond_far_side_sp2"
      #fragment_sets.push *near_side_double_bond_break(carbon, oxygen)
      #end
      #end
      ## note: the case of a carboxy is found with CO search
      #end
      #end
      #end
      #end
      #if rules.include?(:sp3c_nitrogen_double_bond)
      #self.each_match("CN") do |_atoms|
      #(carbon, nitrogen) = _atoms
      #num_nitrogen_bonds = nitrogen.bonds.size
      #case num_nitrogen_bonds
      #when 2
      #if carbon.type == 'C3'
      #fragment_sets.push *near_side_double_bond_break(carbon, nitrogen)
      #end
      #end
      #end
      #end

      #unless had_hydrogens
      #fragment_sets.each {|set| set.each(&:remove_h!) }
      #self.remove_h!
      #end
      #if opts[:uniq]
      ## TODO: impelent properly
      ##fragment_sets = fragment_sets.uniq_by(&:csmiles) 
      #raise NotImplementedError
      #end
      #fragment_sets
      #end

    end
    include Fragmentable
  end
end

