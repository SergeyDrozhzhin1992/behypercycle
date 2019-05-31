clc;

set(0,'DefaultTextInterpreter', 'latex');

set(0,'DefaultAxesFontSize', 24);
set(0,'DefaultTextFontSize', 24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Read custom settings  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_file = fopen('settings_matlab.txt', 'r');
count_type = fread(set_file, 1, 'int');          %The Number of types(sizeA)
count_iter = fread(set_file, 1, 'int') + 1;      %The Number of iterations
count_step = fread(set_file, 1, 'double') + 1;   %The Number of points in the grid to solve ODE
solve_step = fread(set_file, 1, 'int');          %ODU decision step (at which iterations of evolution we will solve an ODE) 
count_solve_step = fread(set_file, 1, 'double'); %The number of iterations, which decided ODE
fclose(set_file);

time = 0 : (count_iter - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read fitness
fitness_file = fopen('fitness_matlab.txt', 'r');
fitness = zeros(count_iter, 1);
for i = 1 : count_iter
    fitness(i, 1) = fread(fitness_file, 1, 'double'); 
end
fclose(fitness_file);

%We draw dependence of fitness on time
figure();
grid on;
hold on;
plot(time, fitness, 'k--', 'LineWidth', 4.5);

xlabel('$\tau$','FontSize', 32);
ylabel('$\bar f$', 'FontSize', 32);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_file = fopen('freqType_matlab.txt', 'r');
freq = zeros(count_type, count_iter);
for i = 1 : count_type
    for j = 1 : count_iter
        freq(i, j) = fread(freq_file, 1, 'double');
    end
end
fclose(freq_file);

% Draw the frequency dependence of time
figure();
grid on;
hold on;

plot(time, freq, 'LineWidth', 3.5);

hold off;

leg = 'u1  ';
for i = 2 : count_type
    if i < 10
        leg = [leg; strcat('u', num2str(i), 32, 32)];
    else if 10 <= i & i < 100
            leg = [leg; strcat('u', num2str(i), 32)];
        else
            leg = [leg; strcat('u', num2str(i))];
        end
    end
end
legend(leg);

xlabel('$\tau$', 'FontSize', 32);
ylabel('$\bar u$', 'FontSize', 32);
%title('Fixed point');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Read the grid in time
time_file = fopen('time_matlab.txt', 'r');
time_cont = zeros(count_step, 1);
for i = 1 : count_step
    time_cont(i, 1) = fread(time_file, 1, 'double');
end    
fclose(time_file);

%Read solutions to the ODE system at every step
freq_cont_file = fopen('freqType_continuos_matlab.txt', 'r');
freq_cell = {};

for i = 1 : count_solve_step
    
    freq = zeros(count_type, count_step);
    for j = 1 : count_type
       for k = 1 : count_step
          freq(j, k) = fread(freq_cont_file, 1, 'double'); 
       end
    end
    
    freq_cell{i} = freq;
end

fclose(freq_cont_file);
    
for i = 1 : count_solve_step
    
    freq = freq_cell{i};  
    
    figure;
    grid on;
    hold on;
    
    plot(time_cont, freq, 'LineWidth', 3);
    
    leg = 'u1  ';
    for j = 2 : count_type
        if j < 10
            leg = [leg; strcat('u', num2str(j), 32, 32)];
        else if 10 <= j & j < 100
                leg = [leg; strcat('u', num2str(j), 32)];
            else
                leg = [leg; strcat('u', num2str(j))];
            end
        end
    end
    legend(leg);

    xlabel('t', 'FontSize', 32);
    ylabel('$\bar u$', 'FontSize', 32);
    
    if i == 1
        title('ODE solution at the 0th iteration');
    else if i == count_solve_step
        title(strcat('ODE solution at the  ', num2str(count_iter - 1), '  iteration'));
        else 
            title(strcat('ODE solution at the  ', num2str((i - 1) * solve_step - 1), '  iteration'));
        end
    end

end

