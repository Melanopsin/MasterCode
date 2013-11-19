function Avg_Matrix_Normalized = Avg_Matrix()

no_runs = 1;
no_species = 9;
dt = 0.01; final_time = floor(15/dt); sizing = 1;

%new_matrx = zeros(size(model_matrx));

new_matrx = zeros(final_time,no_species);
calc_matrx = zeros(no_runs*final_time,no_species);

for j = 1:no_runs
    [t M X] = Melanopsin2(0);
    model_matrx = [t M X];
    
    %model_matrx = sort(10*rand(4));
    %disp('model_matrx = ')

    steps = dt;
    disp(final_time)
    for k = 1:final_time
        new_matrx(k,:) = model_matrx(max(find(model_matrx(:,1)<steps)),:);
        steps = steps + dt;
    end

    calc_matrx(sizing:sizing - 1 + length(new_matrx(:,1)),:) = new_matrx;
    sizing = sizing + length(new_matrx(:,1));
end

t_end = length(new_matrx(:,1));

avging_matrx = zeros(t_end,no_species);

for j = 1:no_runs
    avging_matrx_temp = calc_matrx((j-1)*t_end + 1:j*t_end,:);
    avging_matrx = avging_matrx + avging_matrx_temp;
end

Avg_Matrix = avging_matrx/no_runs;

Normalize_Factor = max(Avg_Matrix);

for i = 2:no_species
    Avg_Matrix(:,i)/Normalize_Factor(i);
end

Avg_Matrix_Normalized = Avg_Matrix;

%plot(Avg_Matrix_Normalized(:,1),Avg_Matrix_Normalized(:,end))

    