function [slopes] = calculate_slopes(time, force, contractpoints, relaxpoints)

        for i = 1:length(contractpoints)
            slopes(i,1:2) = slopeFit(time(contractpoints{i,2}), force(contractpoints{i,2}),0);
            slopes(i,3:4) = slopeFit(time(relaxpoints{i,2}), force(relaxpoints{i,2}),0);
        end
end

function [cfs] =  slopeFit(x,y,phandle)

if length(x) >2
    ok_ = isfinite(x) & isfinite(y);
    if ~all( ok_ )
        warning( 'GenerateMFile:IgnoringNansAndInfs',...
            'Ignoring NaNs and Infs in data.' );
    end
    ft_ = fittype('poly1');
    
    cf_ = fit(x(ok_),y(ok_),ft_);
    
    cfs = coeffvalues(cf_);
   
else
    cfs = NaN;
end

end

