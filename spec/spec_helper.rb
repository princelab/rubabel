require 'rspec'
require 'stringio'

# Requires supporting files with custom matchers and macros, etc,
# in ./support/ and its subdirectories.
#Dir["#{File.dirname(__FILE__)}/support/**/*.rb"].each {|f| require f}

RSpec.configure do |config|
  config.formatter = :documentation  
  config.color = true
end

TESTFILES = File.dirname(__FILE__) + "/testfiles"

module Kernel
  # from: http://thinkingdigitally.com/archive/capturing-output-from-puts-in-ruby/
  def capture_stdout
    out = StringIO.new
    $stdout = out
    yield
    return out.string
  ensure
    $stdout = STDOUT
  end
end
