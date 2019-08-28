function result = extractdata(C,i,type)
% EXTRACTSTR converts the string from the data file into numbers/letters
%   C = importdata(filename)
%   i = integer, index of row of the datafile
%   type = str, 'num' to return a number, 'str' to return a string

    str = split(C{i});
    switch type
        case 'num'
            result = str2double(str(1));
        case 'str'
            result = str(1);
        case 'func'
            str1 = string(str(1));
            strfnc = ['@(x,y) ', str1];
            joined = join(strfnc);
            result = str2func(joined);
    end
end