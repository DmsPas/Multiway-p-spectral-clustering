function [Cut,CCut] = computeRCutValue(clusters,W,normalized)
% Computes the components in the Ratio/Normalized Cut and Ratio/Normalized Cheeger Cut expression.
%
% Usage: [cutpart1,cutpart2] = computeCutValue(clusters,W,normalized)

K = max(clusters);

% Initialized Cut and Cheeger Cut
Cut  = 0;
CCut = 0;

for k = 1:K
    
    W2  = W(clusters==k,clusters~=k);
    cut = full(sum(sum(W2)));
    
    
    if (~normalized)
        cardinalityA = sum(clusters==k);
        cardinalityB = sum(clusters~=k);
        
        cutpart = cut/cardinalityA;
        cutpart_min = cut/min(cardinalityA,cardinalityB);
    else
        degreeA = sum(W(:,clusters==k));
        degreeB = sum(W(:,clusters~=k));
        
        volA   = sum(degreeA);
        volB   = sum(degreeB);
        
        cutpart = cut/volA;
        cutpart_min = cut/min(volA,volB);       
    end
    Cut  = Cut  + cutpart; 
    CCut = CCut + cutpart_min;  
end

end
