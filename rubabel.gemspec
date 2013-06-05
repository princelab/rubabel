# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'rubabel/version'

Gem::Specification.new do |spec|
  spec.name          = "rubabel"
  spec.version       = Rubabel::VERSION
  spec.authors       = ["John Prince"]
  spec.email         = ["jtprince@gmail.com"]
  spec.description   = %q{TODO: Write a gem description}
  spec.summary       = %q{TODO: Write a gem summary}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files`.split($/)
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  [
    ["openbabel", "~> 2.3.1.2"],
    ["andand", "~> 1.3.3"],
    ["mini_magick", "~> 3.6.0"]
  ].each do |args|
    gem.add_dependency(*args)
  end
  [
    ["bundler", "~> 1.3"],
    ["rspec", "~> 2.13.0"], 
    ["rdoc", "~> 3.12"], 
    ["jeweler", "~> 1.8.3"],
    ["rake"]
  ].each do |args|
    spec.add_development_dependency *args
  end
end
