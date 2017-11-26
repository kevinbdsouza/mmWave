%getCDLProfile TR 38.900 V14.2.0 Section 7.7.1 Clustered Delay Line (CDL) per-cluster parameters

%   Copyright 2016-2017 The MathWorks, Inc.

function out = getCDLPerClusterParameters(delayProfile)

    persistent CDL_A CDL_B CDL_C CDL_D CDL_E;
    uXprUmaNlos = 7;
    uXprUmaLos = 8;
        
    if (isempty(CDL_A))

        % Table 7.7.1-1
        CDL_A = makePerClusterParameters( 5, 11,  3,  3, uXprUmaNlos);
                
        % Table 7.7.1-2
        CDL_B = makePerClusterParameters(10, 22,  3,  7,  8);
        
        % Table 7.7.1-3
        CDL_C = makePerClusterParameters( 2, 15,  3,  7,  7);
        
        % Table 7.7.1-4
        CDL_D = makePerClusterParameters( 5,  8,  3,  3, uXprUmaLos);
        
        % Table 7.7.1-5
        CDL_E = makePerClusterParameters( 5, 11,  3,  7,  8);

    end

    switch (upper(delayProfile))

        case 'CDL-A'
            out = CDL_A;
            %disp('hello');
        case 'CDL-B'
            out = CDL_B;
        case 'CDL-C'
            out = CDL_C;
        case 'CDL-D'
            out = CDL_D;
        case 'CDL-E'
            out = CDL_E;

    end

    function out = makePerClusterParameters(C_ASD,C_ASA,C_ZSD,C_ZSA,XPR)

            out.C_ASD = C_ASD;
            out.C_ASA = C_ASA;
            out.C_ZSD = C_ZSD;
            out.C_ZSA = C_ZSA;
            out.XPR = XPR;

    end
    
end
