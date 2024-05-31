function s = csv2struct(filename)

   % Get the header
   fid = fopen(filename, 'r');
   if fid < 0
      error(['Can''t open file "' filename '" for reading.']);
   end
   header = fgetl(fid); % Read header line
   fclose(fid);

   % Convert header string to a cell array with field names.
   fields = eval(['{''', strrep(header, ',', ''','''), '''}']);

   % Get the data
   data = csvread(filename, 1, 0);

   % Convert data into a cell array of values.
   values = num2cell(data, 1);

   % Build structure
   list = { fields{:} ; values{:} };
   s = struct(list{:});