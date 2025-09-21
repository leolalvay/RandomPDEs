% I will plot the old density and the new density in 3D and the contour
% plots
%next try to see what I did in the code and what I can add from it

%As far as I can remember, h is the function iside the integral (rho*f) and
%the goal of this is to see if its integral is truely 0 (check the
%crediblity of crude MC to estimate I)

K_values = [3, 6, 12, 18, 24, 30];
num_values = numel(K_values);

% Create a figure
figure;

% Iterate over the K values and plot h in subplots
for i = 1:num_values
    K = K_values(i);
    
    % Create a subplot for each K value
    subplot(2, 3, i); % Adjust the grid size according to your preference
    
    % Call the function plot_h(K) to plot h
    plot_h(K);
end







%%
%%%%%%FUNCTIONS%%%%%%


%% The function f
function result = f(x, K)
    result = max(exp(x(:, 1)) + exp(x(:, 2)) - K, 0);
end

%%The function h inside of the integral
function result= h(x,K) %the function rho.f to find u1* and u2*
        % result= f(x,K) .* normpdf(x(:,1), 0, 1) .* normpdf(x(:,2), 0, 1);
        result= f(x,K) .*rho(x);
end

%% The density rho (standard normal of x=(x1,x2)
 function result=rho(x)
        result= normpdf(x(:,1), 0, 1) .* normpdf(x(:,2), 0, 1); %!(1/sqrt(2*pi)) is already inside the normpdf (in the expression of the normal pdf) and you should not add it
 end

%% The function to plot h for a value of K
% Function to plot h in 3D
function plot_h(K)
    x_ax1 = linspace(-5, 5, 100);
    x_ax2 = linspace(-5, 5, 100);
    [X_ax1, X_ax2] = meshgrid(x_ax1, x_ax2);
    x = [X_ax1(:), X_ax2(:)];
    
    % Calculate the values of h
    H = h(x, K);
    H = reshape(H, size(X_ax1));
    
    % Plot the 3D function
    surf(X_ax1, X_ax2, H);
    colorbar;
    xlabel('X1');
    ylabel('X2');
    zlabel('h');
    titleString = sprintf('h for K = %d', K);
    title(titleString);
end