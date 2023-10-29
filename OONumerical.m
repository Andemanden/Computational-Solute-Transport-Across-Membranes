%%objectoriented numerical method
classdef OONumerical
    properties
        domain_length
        domain_steps
        dx
        time_length
        time_steps
        dt
        D
        velocity
        feed_conc
        rejection_rate
        x
        t
        C
    end
    
    methods
        function obj = OONumerical(domain_length, domain_steps, time_length, time_steps, D, velocity, feed_conc, rejection_rate)
            obj.domain_length = domain_length;
            obj.domain_steps = domain_steps;
            obj.dx = domain_length / (domain_steps - 1);
            obj.time_length = time_length;
            obj.time_steps = time_steps;
            obj.dt = time_length / time_steps;
            obj.D = D;
            obj.velocity = velocity;
            obj.feed_conc = feed_conc;
            obj.rejection_rate = rejection_rate;
            
            % Create a grid for space and time
            obj.x = linspace(0, domain_length, domain_steps);
            obj.t = linspace(0, time_length, time_steps);
            
            % Initialize the concentration array (1D)
            obj.C = ones(domain_steps, time_steps);
        end
        
        function simulate(obj)
            for j = 2:obj.time_steps
                for i = 1:obj.domain_steps
                    if i == 1
                        obj.C(1, :) = obj.feed_conc;
                        d2Cdx2 = obj.D * (obj.C(2, j-1) - obj.C(1, j-1)) / obj.dx^2;
                        dCdt_advection = -obj.velocity * obj.C(i, j-1) / obj.dx;
                    elseif i == obj.domain_steps
                        d2Cdx2 = (obj.D * (obj.C(obj.domain_steps - 1, j-1) - obj.C(obj.domain_steps, j-1)) / obj.dx^2) + obj.C(i, j-1) * (obj.rejection_rate) * obj.velocity / (obj.dx);
                    else
                        d2Cdx2 = obj.D * (obj.C(i + 1, j-1) - 2 * obj.C(i, j-1) + obj.C(i - 1, j-1)) / obj.dx^2;
                        dCdt_advection = -obj.velocity * (obj.C(i, j-1) - obj.C(i-1, j-1)) / obj.dx;
                    end
                    obj.C(i, j) = obj.C(i, j-1) + obj.dt * (d2Cdx2 + dCdt_advection);
                end
            end
        end
        
        function plotResults(obj, time_fraction)
            % Define fractions of time steps you want to visualize
            time_instances = round(time_fraction * obj.time_steps);
            
            % Create a figure for the first plot
            figure;
            hold on;
            
            % Plot concentration at specific time instances
            for i = 1:length(time_instances)
                time_index = time_instances(i);
                plot(obj.x, obj.C(:, time_index));
            end
            
            xlabel('Position (meters)');
            ylabel('Concentration');
            title('Concentration Over Time at Different Instances');
            
            % Add a legend for clarity
            legend(arrayfun(@(f) ['t=', num2str(f)], time_fraction, 'UniformOutput', false));
            
            grid on;
            hold off;
            
            % Create a 3D surface plot to visualize concentration over time and position
            [T, X] = meshgrid(obj.t, obj.x);
            figure;
            h = surf(X, T, obj.C);
            xlabel('Position (meters)');
            ylabel('Time (seconds)');
            zlabel('Concentration');
            title('Concentration Over Time and Position');
            
            % Set axis limits to start at the origin
            xlim([0, obj.domain_length]);
            ylim([0, obj.time_length]);
            zlim([0, max(obj.C(:))]);
            set(h, 'LineStyle', 'none');
        end
    end
end