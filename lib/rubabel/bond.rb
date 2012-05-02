
module Rubabel
  class Bond

    attr_accessor :obbond

    def initialize(obbond)
      @obbond
    end

    # does it include the Rubabel::Bond
    def include?(atom)
    end

    # does it include the OpenBabel::OBAtom
    def ob_include?(obatom)
    end

    %w(aromatic in_ring rotor amide primary_amide secondary_amide tertiary_amide
       ester carbonyl single double triple ksingle kdouble ktriple closure
       up down wedge hash wedge_or_hash cis_or_trans double_bond_geometry
      ).each do |bool|
      # can return true, false or nil
      define_method("#{bool}?") do
        @obbond.send("is_#{bool}")
      end

      # can be set to nil, true or false
      define_method("#{bool}?=") do |val|
        @obbond.send("set_#{bool}", val)
      end
    end
  end
end
