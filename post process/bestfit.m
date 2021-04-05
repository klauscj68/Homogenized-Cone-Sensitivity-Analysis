function [varargout] = bestfit(x,y,flag_model)
%Find the best fitting quadratic or exponential to the data
%   x is the vector of indvar samples. y is the vector of depvar samples.
%   flag_model can be either 'quad' or 'exp' to indicate what we use to
%   fit. varargout is the matlab call for outputting a possibly variable
%   number of of arguments.  Basically, nargout is the number of outputs
%   the user told the code to output when they entered the line [...] =
%   bestfit(...), and that's how matlab knows how many arguments to return.
%   NOTE: Code uses the curve fitting toolbox to avoid writing custom
%         Newton-Raphson code

switch flag_model
    case 'quad'
        % Check if right number of outputs
        if nargout ~= 1
            error('For quadratic model, there is one output coeff')
        end
        
        % Define the .5*A*x^2 model.
        %  Specify how the best fit will be performed using fitoptions. Do
        %  this by least squares which is known to be linear in the unknown
        %  coefficient, ie adjustable parameter, A
        fopt = fitoptions('METHOD','LinearLeastSquares');
        
        %  Since it is linear in A, we can use a cell-based EXPR to encode
        %  the fittype.  
        ft = fittype({'.5*x^(2)'},'coefficients','A','options',fopt);
        
        %  Perform the fit
        fobj = fit(x(:),y(:),ft);
        
        varargout{1} = fobj.A;
        
    case 'exp'
        % Check if right number of outputs
        if nargout ~= 2
            error('For exponential model, there are two output coeffs')
        end
        
        % Check if y-values are all positive
        check = sum(y <= 0);
        if check ~= 0
            disp('You are trying to fit an exponential to nonpositives. Performing shift.')
            
            % Shift data by least value
            y = y - (1.0001)*min(y);
        end
        
        % Define and fit the c*exp(-alpha*t) model. Do this by taking a ln
        %  and then fitting a line.
        y = log(y);
        
        fopt = fitoptions('METHOD','LinearLeastSquares');
        ft  = fittype({'1','x'},'coefficients',{'ln_c','alpha'},...
                      'options',fopt);
        
        fobj = fit(x(:),y(:),ft);
        varargout{1} = fobj.alpha;
        varargout{2} = exp(fobj.ln_c);
        
    otherwise
        error('flag_model must either be quad or exp')
end

end

