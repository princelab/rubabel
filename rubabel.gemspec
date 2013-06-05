# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'rubabel/version'

Gem::Specification.new do |spec|
  spec.name          = "rubabel"
  spec.version       = Rubabel::VERSION
  spec.authors       = ["John T. Prince"]
  spec.email         = ["jtprince@gmail.com"]
  spec.description   = %q{Ruby interface to the OpenBabel ruby bindings similar to pybel}
  spec.summary       = %q{Ruby interface to the openbabel ruby bindings (or the openbabel gem).  The
interface attempts to be a ruby-ish analogue of pybel.}
  spec.homepage      = "http://github.com/princelab/rubabel"
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
    spec.add_dependency(*args)
  end
  [
    ["bundler", "~> 1.3"],
    ["rspec", "~> 2.13.0"], 
    ["rdoc", "~> 3.12"], 
    ["jeweler", "~> 1.8.3"],
    ["simplecov"],
    ["rake"]
  ].each do |args|
    spec.add_development_dependency(*args)
  end
end
