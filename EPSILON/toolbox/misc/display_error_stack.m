function [] = display_error_stack(err)
    disp('')
    for ii=1:length(err.stack)
        disp(['error in ' err.stack(ii).name ' - line ' num2str(err.stack(ii).line)])
    end
    error('')