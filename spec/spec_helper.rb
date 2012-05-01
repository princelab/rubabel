require 'rspec'

# Requires supporting files with custom matchers and macros, etc,
# in ./support/ and its subdirectories.
#Dir["#{File.dirname(__FILE__)}/support/**/*.rb"].each {|f| require f}

RSpec.configure do |config|
  config.formatter = :documentation  
  config.color = true
end

TESTFILES = File.dirname(__FILE__) + "/testfiles"
