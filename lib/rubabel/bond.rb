
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

  end
end
