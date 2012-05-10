
module Rubabel
  class Fragmentation

    ADDUCTS = [:lioh, :nh4cl, :nh4oh]

    DEFAULTS = {
      ph: 7.4,
      adducts: [],
    }

    def initialize(mol, options)
      @mol = mol
      @options = options
    end

    def fragment(rules=[:co])
    end



  end
end
