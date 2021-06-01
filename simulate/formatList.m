function fl = formatList(l, f)
    cl = arrayfun(@(x) sprintf(f, x), l, 'UniformOutput', false);
    fl = strjoin(cl, ', ');
end
