%% Rename field in a structure

function newfield = renamefield(mystruct, newnames)
    F = fieldnames(mystruct);
    newfield = struct;
    for i = 1:length(F)
        newfield.(newnames{i}) = mystruct.(F{i});
    end
end