# encoding: utf-8

require 'rubygems'
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "rubabel"
  gem.homepage = "http://github.com/princelab/rubabel"
  gem.license = "MIT"
  gem.summary = %Q{Ruby interface to the OpenBabel ruby bindings similar to pybel}
  gem.description = %Q{Ruby interface to the openbabel ruby bindings (or the openbabel gem).  The
interface attempts to be a ruby-ish analogue of pybel.}
  gem.email = "jtprince@gmail.com"
  gem.authors = ["John T. Prince"]
  [
    ["shoulda", ">= 0"], 
    ["rdoc", "~> 3.12"], 
    ["jeweler", "~> 1.8.3"]
  ].each do |gem, version_string|
    gem.add_development_dependency(gem, version_string)
  end
end
Jeweler::RubygemsDotOrgTasks.new

require 'rspec/core'
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

task :default => :spec

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "rubabel #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
