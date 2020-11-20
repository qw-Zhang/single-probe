function PAS = generate_PAS(ideal_phi,pow,AS,step_mov)
%ideal scenario
amp = 10.^(pow/10);
A = amp./max(amp);
PAS_norm = zeros(1,length(ideal_phi));
for i = 1:length(ideal_phi)
    % Lap distr
    PAS_norm(i) = 1*(1/(2*AS))*(exp(-abs(ideal_phi(i)-0)/AS));
end

PAS_move_temp = zeros(length(step_mov),length(PAS_norm));
for j = 1:length(step_mov)
    PAS_norm = A(j)*PAS_norm;
    if(step_mov(j) < 0)
        PAS_move_temp(j,:) = [PAS_norm(abs(step_mov(j))+1:length(PAS_norm)),PAS_norm(1:abs(step_mov(j)))];
    else
        PAS_move_temp(j,:) = [PAS_norm(length(PAS_norm) - abs(step_mov(j)) + 1:length(PAS_norm))...
            ,PAS_norm(1:length(PAS_norm) - abs(step_mov(j)))];
    end
end
PAS = sum(PAS_move_temp,1);
%using for normalized PAS
sum_PAS = sum(PAS);
PAS = PAS/sum_PAS;