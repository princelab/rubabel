
# puts if $VERBOSE
def putsv(*args)
  if $VERBOSE
    puts(*args)
  end
end

# puts to $stderr if $VERBOSE
def putsev(*args)
  if $VERBOSE
    $stderr.puts(*args)
  end
end
