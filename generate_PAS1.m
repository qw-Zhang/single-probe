% generate PAS data by function:PAS_fun
function PAS = generate_PAS1(ideal_phi,scenario)
    switch scenario
        case 'test'
            AOA = -22.5*pi/180;AS = 35*pi/180;
            pow = 0;
            
        case 'micro'
            % % SCME urban micro-cell
            AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;
            pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
            
        case 'macro'
            % %SCME urban macro-cell
            AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;AS = 35*pi/180;
            pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
            
        otherwise
            warning('wrong scenario.');
            quit;
    end
    
    amp = 10.^(pow/10);
    A = amp./max(amp);
    PAS_fun = @(x) (1/(2*AS))*sum( A.*exp(-abs(x-AOA)/AS) );
    
    PAS = zeros(1,length(ideal_phi));
    for i = 1:length(ideal_phi)
        PAS(i) = PAS_fun(ideal_phi(i));
    end
    sumPAS = sum(PAS);
    PAS = PAS./sumPAS;
end