#!/usr/bin/env ruby

require 'trollop'
require 'rubabel'
require 'rubabel/molecule/fragmentable'


parser = Trollop::Parser.new do
  banner "usage: #{File.basename($0)} [OPTIONS|RULES] <SMARTS> ..."
  text "\noptions:"
  opt :ph, "the pH to use (experimental option)", :default => Rubabel::Molecule::Fragmentable::DEFAULT_OPTIONS[:ph]
  #opt :uniq, "no repeated fragments", :default => false
  text "\nrules:"
  Rubabel::Molecule::Fragmentable::RULES.each do |rule|
    opt rule, rule.to_s.gsub("_",' ')
  end
  text "\nexample:"
  text "fragmenter.rb -aecsoxn 'CCC(=O)OCCC' 'CCC(=O)OCCC(=O)O'"
end

rules = parser.parse(ARGV)
options = {rules: []}
options[:ph] = rules.delete(:ph)
options[:uniq] = rules.delete(:uniq)
rules.each do |k,v|
  options[:rules] << k if v && k.to_s !~ /_given/
end

if ARGV.size == 0
  parser.educate && exit
end

ARGV.each do |mol|
  mol = Rubabel[mol]
  puts "\nmolecule: #{mol.csmiles}"
  fragment_sets = mol.fragment(options)
  fragment_sets.each do |frag_set|
    puts ""
    frag_set.each do |frag|
      puts "#{frag.mass.round(5)} #{frag.csmiles}"
    end
  end
end
