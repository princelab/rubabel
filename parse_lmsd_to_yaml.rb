require 'yaml'
require 'pry'

file = 'LMSD_20120412_All_sd.csv'

lines = File.readlines(file)
header = lines.shift.delete('"').split(",")

keys = [:pubchem_url, :lmaps_url, :lmid, :systematic_name, :synonyms, :lipid_category, :lipid_class, :exact_mass, :formula, :pubchem_sid, :inchi_key, :common_name, :kegg_id, :chebi_id, :subclass, :hmdbid, :lipidbank_id, :lipidat_id, :metabolomics_id, :class_level4, :structure]

LMID = {}

lines.map do |line|
  arr = line.split('","')
  begin
    row = Hash[keys[2] => arr[2], keys[-1] => arr[-1]]
  raise ArgumentError if row[:structure].nil?
  rescue ArgumentError => e
    binding.pry
    p e
    p arr
    p arr.size
    abort
  end
  row[:structure] = row[:structure].split("|").join("\n")
  LMID[row[:lmid]] = row
end

# Write all LMIDS to file
File.open("all_lmids.txt",'w') do |io|
  io.write LMID.keys
end

File.open("sdfs.yml",'w') do |io|
  io.write LMID.to_yaml
end
