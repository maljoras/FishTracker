function opts = setfield(opts,field,value)

  eval(sprintf('opts.%s = %f;',field,value));