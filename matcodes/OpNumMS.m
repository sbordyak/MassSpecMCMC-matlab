function n=OpNumMS(oper)

if strcmp(oper,'changer')
    n=1;
elseif strcmp(oper,'changeI')
    n=2;
elseif strcmp(oper,'changedfg')
    n=3;
elseif strcmp(oper,'changebl')
    n = 4;
elseif strcmp(oper,'noise')
    n = 5;
end