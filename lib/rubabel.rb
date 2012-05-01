require 'openbabel'

%w(atom molecule fingerprint smarts molecule_data).each do |klass|
  require "rubabel/#{klass}"
end

module Rubabel
  # determines the extension from filename if type is nil
  def self.foreach(filename, type=nil, &block)
    block or return enum_for(__method__, filename, type)
    type ||= File.extname(filename)[1..-1]
    obconv = OpenBabel::OBConversion.new
    obconv.set_in_format(type.to_s)
    obmol = OpenBabel::OBMol.new
    # gets the first mol
    not_at_end = obconv.read_file(obmol, filename)
    while !not_at_end
      block.call Rubabel::Molecule.new(obmol)
      obmol = OpenBabel::OBMol.new
      not_at_end = obconv.read_file(obmol, filename)
      p filename
      p not_at_end
    end
  end
end


=begin
OBConversion conv;
OBMol mol;
bool success = conv.SetInFormat("sdf");
if(success)
  {
    bool notatend = conv.ReadFile(&mol, "myfile.sdf");
    // Do something with mol
    while(notatend)
      {
        notatend = conv.Read(&mol);
        // Do something with mol
      }
  }
=end
