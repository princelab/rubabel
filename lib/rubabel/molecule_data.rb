require 'openbabel'

module Rubabel
  # A hash-like interface for dealing with tag data in molecules inspired by
  # the MoleculeData implementation in pybel.  This does not actually hash the
  # data, it just reads it from the OBMol object each time and acts like a
  # hash.
  class MoleculeData
    include Enumerable

    OB_PAIR_DATA_TYPE = OpenBabel::OBPairData.new.get_data_type
    OB_COMMENT_DATA_TYPE = OpenBabel::OBCommentData.new.get_data_type

    DATA_TYPES = [OB_PAIR_DATA_TYPE, OB_COMMENT_DATA_TYPE]

    attr_accessor :obmol

    def initialize(obmol)
      @obmol = obmol
    end

    def pair_data(&block)
      block or return enum_for(__method__)
      # TODO: should do this with an ob iterator
      @obmol.get_data.each do |ob_gd|
        if DATA_TYPES.include?( ob_gd.get_data_type )
          block.call( OpenBabel.to_pair_data( ob_gd ) )
        end
      end
    end

    def each(&block)
      block or return enum_for(__method__)
      pair_data.each do |pd|
        block.call [pd.get_attribute, pd.get_value]
      end
    end

    def to_a
      each.to_a
    end

    def size
      pair_data.to_a.size
    end
    alias_method :length, :size

    def [](key)
      pair_data.find {|pd| pd.get_attribute == key }.get_value
    end

    # returns the val
    def []=(key,val)
      if key?(key)
        OpenBabel.to_pair_data(@obmol.get_data(key)).set_value(val)
        val
      else
        pd = OpenBabel::OBPairData.new
        pd.set_attribute(key)
        pd.set_value(val)
        @obmol.clone_data(pd)
        val
      end
    end

    def key?(key)
      pair_data.any? {|pd| pd.get_attribute == key }
    end

    def keys
      pair_data.map(&:get_attribute)
    end

    def values
      pair_data.map(&:get_value)
    end

    def delete(key, &block)
      if key?(key)
        val = self[key]
        @obmol.delete_data( @obmol.get_data(key) )
        val
      else
        block ? block.call : nil
      end
    end

  end
end
