function result = extractstr(C,i,type)
% EXTRACTSTR converts the string from the data file into numbers/letters
%   C = importdata(filename)
%   i = integer, index of row of the datafile
%   type = str, 'num' to return a number, 'str' to return a string

    str = split(C{i});
    if type == 'num'
        result = str2double(str(1));
    elseif type == 'str'
        result = str(1);
    end
end
