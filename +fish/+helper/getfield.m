function value = getfield(opts,field)

  eval(sprintf('value = opts.%s;',field));