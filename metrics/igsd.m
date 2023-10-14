
myFunc = @(u, v) (0.25*sin(4*u)*sin(v)*(cos(-4*v) + cos(-v))) - 0.25*(sin(2*u)*sin(2*u)*(cos(-5*v)-cos(-3*v))) - (sin(4*u) * sin(-4*v) *(cos(v)+1) ) - (0.5*sin(4*u) *sin(-v));

xmin = -1;
xmax = 1;
ymin = -1.5;
ymax = -1;
num_samples = 200;
h = 0.001; % Adjust the step size as needed

% Call the function to compute the standard deviation of the derivative
std_dev_derivative = computeStdDeviationOfDerivative(myFunc, xmin, xmax, ymin, ymax, num_samples, h);

% Display the result
fprintf('IGSD: %.4f\n',1/std_dev_derivative);

function std_deviation_derivative = computeStdDeviationOfDerivative(func, xmin, xmax, ymin, ymax, num_samples, h)
    % Input arguments:
    % func: A function handle representing the two-dimensional function.
    % xmin, xmax, ymin, ymax: The range for x and y values for sampling.
    % num_samples: The number of samples to generate.
    % h: Step size for numerical differentiation.

    % Generate random samples within the specified range
    x_samples = xmin + (xmax - xmin) * rand(1, num_samples);
    y_samples = ymin + (ymax - ymin) * rand(1, num_samples);

    % Initialize an array to store the derivative values at the sampled points
    derivative_samples = zeros(1, num_samples);

    % Evaluate the derivative of the function at the sampled points
    for i = 1:num_samples
        x = x_samples(i);
        y = y_samples(i);
        
        % Use central difference for numerical derivative
        func(x, y+h)
        func(x,y-h)
        derivative_samples(i) = (func(x+h, y) - func(x-h,y)) / (2 * h);
    end
    
    % Compute the standard deviation of the derivative values
    std_deviation_derivative = std(derivative_samples);
end
