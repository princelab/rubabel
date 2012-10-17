#!/usr/bin/env ruby

require 'trollop'
require 'rubabel'
require 'rubabel/molecule/fragmentable'

default_ph = 2.5

parser = Trollop::Parser.new do
  banner "usage: #{File.basename($0)} [OPTIONS|RULES] <SMARTS> ..."
  text "\noptions:"
  opt :ph, "the pH to use (experimental option)", :default => default_ph
  opt :images, "print out svg images of fragments" 
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
  mol = Rubabel[smiles]
  puts "\nmolecule: #{mol.csmiles}"
  mol.correct_for_ph!(ph)
  puts "at ph #{ph}: #{mol.csmiles}"
  fragment_sets = mol.fragment(opts)
  puts %w(mass charge smiles pairing).join("\t")
  fragment_sets.each_with_index do |frag_set,i|
    frag_set.each_with_index do |frag,j|
      if options[:images]
        frag.title = 
          if frag.charge == 0
            'nocharge'
          else
            mz=(frag.mass / frag.charge).round(5).to_s
          end
        fn = "frag#{i}-#{j}_#{frag.title}.png"
        frag.write(fn)
      end
      puts [frag.mass.round(5), frag.charge, frag.csmiles, i].join("\t")
    end
  end
end
