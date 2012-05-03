module ObjectSpace
  class << self
    include Enumerable
    alias_method :each, :each_object
  end
  def self.modules
    self.select {|obj| obj.is_a?(Module) }
  end
end

class Dog
end

before = ObjectSpace.modules

require 'openbabel'

after = ObjectSpace.modules

modules = (after - before).reject {|mod| mod.to_s =~ /SWIG/i }

class Object
  def arity_hash
    uniq = self.methods - (Object.new.methods + Object.methods)
    Hash[ uniq.map do |name|
      [name.to_s, self.method(name).arity]
    end
    ]
  end
end

modules.each do |mod|
  obj_hash = 
    if mod.respond_to?(:allocate) 
      begin
        obj = mod.allocate
        obj.arity_hash
      rescue
        {}
      end
    else
      {}
    end
  klss_hash = mod.arity_hash
  if (obj_hash.size + klss_hash.size) > 0
    puts "*" * 50
    puts mod
    puts "*" * 50
    puts obj_hash.map {|pair| pair.join(" ") }
    puts klss_hash.map {|pair| pair.join(" ") }
  end
end
