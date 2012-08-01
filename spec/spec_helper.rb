require 'rspec'
require 'stringio'

# Requires supporting files with custom matchers and macros, etc,
# in ./support/ and its subdirectories.
#Dir["#{File.dirname(__FILE__)}/support/**/*.rb"].each {|f| require f}

RSpec.configure do |config|
  config.treat_symbols_as_metadata_keys_with_true_values = true
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

  def capture_stderr
    out = StringIO.new
    $stderr = out
    yield
    return out.string
  ensure
    $stderr = STDERR
  end


end
