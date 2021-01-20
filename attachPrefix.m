function new_labels = attachPrefix(prefix, labels)
% ATTACHPREFIX
%
%   new_labels = attachPrefix(prefix, labels)
%
%   Assuming a cell array of strings, labels, we attach the string,
%   prefix, to the beginning of each item.
    
    new_labels = cell(size(labels));
    
    for i=1:length(labels)
        new_labels{i} = [prefix labels{i}];
    end

end