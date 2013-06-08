#!/usr/bin/env ruby

require 'trollop'
require 'rubabel'
require 'rubabel/molecule/fragmentable'

default_ph = 2.5

Fragment = Struct.new(:frag, :id, :title, :mz, :mass, :charge, :smiles, :pairing)

parser = Trollop::Parser.new do
  banner "usage: #{File.basename($0)} [OPTIONS|RULES] <SMARTS> ..."
  text "\noptions:"
  opt :ph, "the pH to use (experimental option)", :default => default_ph
  opt :images, "print out svg images of fragments" 
  opt :format, "format of the molecules", :default => 'smiles'
  #opt :uniq, "no repeated fragments", :default => false
  text "\nrules:"
  Rubabel::Molecule::Fragmentable::RULES.each do |rule|
    opt rule, rule.to_s.gsub("_",' ')
  end
  text "\nexample:"
  text "fragmenter.rb -xeh 'CCC(=O)OCCC' 'CCC(=O)OCCC(=O)O'"
end

options = parser.parse(ARGV)
opts = {rules: []}
opts[:uniq] = options.delete(:uniq)
ph = options.delete(:ph)
opts[:rules] = Rubabel::Molecule::Fragmentable::RULES.map do |rule|
  rule if options["#{rule}_given".to_sym]
end.compact

if ARGV.size == 0
  parser.educate && exit
end

ARGV.each do |smiles|
  mol = Rubabel[smiles, options[:format].to_sym]
  puts "\nmolecule: #{mol.csmiles}"
  mol.correct_for_ph!(ph)
  puts "at ph #{ph}: #{mol.csmiles}"
  fragment_sets = mol.fragment(opts)
  puts %w(mz mass charge title smiles pairing).join("\t")
  frags = []
  fragment_sets.each_with_index do |frag_set,i|
    frag_set.each_with_index do |frag,j|
      unless frag.charge == 0
        mz = (frag.mass / frag.charge).round(5)
      end

      frag.title = "#{i}-#{j}pair_" + (mz ? "#{mz}_mz" : "#{frag.mass.round(3)}_Mass")
      frag_obj = Fragment.new(frag, frag.title, frag.title, mz, frag.exact_mass, frag.charge, frag.csmiles, i)
      frags << frag_obj
    end
  end
  frags = frags.sort_by {|frag| [-frag.charge, frag.mz] }
  if options[:images]
    frags.each do |frag|
      fn = "#{frag.title}.svg"
      frag.frag.write(fn)
    end
  end
  frags.each do |frag|
    puts [:mz, :mass, :charge, :title, :smiles, :pairing].map {|cat| frag.send(cat) }.join("\t")
  end
end
