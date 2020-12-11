% generate PAS data by function:PAS_fun
function PAS = generate_PAS1(PAS_fun,ideal_phi)
    PAS = zeros(1,length(ideal_phi));
    for i = 1:length(ideal_phi)
        PAS(i) = PAS_fun(ideal_phi(i));
    end
sumPAS = sum(PAS);
PAS = PAS./sumPAS;