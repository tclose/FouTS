function fig = plot_fibres(varargin)

  if nargin == 0
    error('No arguments supplied');
  end

  filename = varargin{1};

  if file_extension(filename) == 'str'
    fig = plot_strands(varargin{:});
  elseif file_extension(filename) == 'sst'
    fig = plot_strand_sets(varargin{:});
  elseif file_extension(filename) == 'tct'
    fig = plot_tracts(varargin{:});
  elseif file_extension(filename) == 'tst'
    fig = plot_tract_sets(varargin{:});
  elseif file_extension(filename) == 'tck'
    fig = plot_tcks(varargin{:});
  elseif file_extension(filename) == 'kst'
    fig = plot_tck_sets(varargin{:});    
  else
    error(['Unrecognised extension, ''' file_extension(filename) '''.']);
  end
  
end