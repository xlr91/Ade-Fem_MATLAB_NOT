function result = extractstr(C,i,type)
    str = split(C{i});
    if type == 'num'
        result = str2double(str(1));
    elseif type == 'str'
        result = str(1);
    end
end
