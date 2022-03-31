function oper = RandomOperMS(hier)

% Randomly generate next model operation, with or without hierarchical step

if hier == 0
        
        nxt = randi(100);
        
        if nxt<41 
            oper = 'changeI';
        elseif nxt<61            
            oper = 'changer';
        elseif nxt<81
            oper = 'changebl';
        elseif nxt<101
            oper = 'changedfg';
        end
        
        
   
    
elseif hier==1
    
        
        nxt = randi(120);
        
        
        if nxt<61 
            oper = 'changeI';
        elseif nxt<81            
            oper = 'changer';
        elseif nxt<96
            oper = 'changebl';
        elseif nxt<101
            oper = 'changedfg';
        elseif nxt<121
            oper = 'noise';
        end
        
  
    
    
else
    display('Is it hierarchical or not?')
end


