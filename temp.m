
clear all

no_runs = 5;
no_species = 4;
final_time = 3; dt = 5; sizing = 1;
%t_end = final_time/dt;

%new_matrx = zeros(size(model_matrx));
new_matrx = zeros(4);
calc_matrx = zeros(no_runs*length(new_matrx(:,1)),no_species);

for j = 1:no_runs
    %model_matrx = Melanopsin();
    model_matrx = sort(10*rand(4));
    disp('model_matrx = ') 
    model_matrx

    %model_matrx
    steps = dt;
    for k = 1:final_time
        new_matrx(k,:) = model_matrx(max(find(model_matrx(:,1)<steps)),:);
        steps = steps + dt;
    end
disp('new_matrx = ') 
new_matrx
%keyboard;
    calc_matrx(sizing:sizing - 1 + length(new_matrx(:,1)),:) = new_matrx;
    disp('calc_matrx =')
    calc_matrx
    sizing = sizing + length(new_matrx(:,1));
end

t_end = length(new_matrx(:,1));

avging_matrx = zeros(t_end,no_species);

disp('calc_matrx = ')
calc_matrx

for j = 1:no_runs
    avging_matrx_temp = calc_matrx((j-1)*t_end + 1:j*t_end,:);
    avging_matrx = avging_matrx + avging_matrx_temp;
end

disp('avging_matrx = ')
avging_matrx

mean_matrx = avging_matrx/no_runs
    