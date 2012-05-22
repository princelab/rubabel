
module Rubabel
  class Molecule
    module Fragmentable
      RULES = [:co]
      ADDUCTS = [:lioh, :nh4cl, :nh4oh]

      DEFAULT_OPTIONS = {
        rules: RULES,
        adduct: nil,
        ph: 7.4,
      }

      # will add_h!(ph) at the given ph
      def fragment(opts={})
        opts = DEFAULT_OPTIONS.merge(opts)
      end
    end

    class Molecule
      include Fragmentable
    end
  end
end
