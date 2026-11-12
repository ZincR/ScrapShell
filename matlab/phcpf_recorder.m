function varargout = phcpf_recorder(cmd, varargin)
%PHCPF_RECORDER Store the iterates of PH-CPF's Levenberg-Marquardt solve for the viewer.
%   phcpf_recorder('start')              begin recording (clears previous data)
%   phcpf_recorder('stop')               stop recording
%   tf = phcpf_recorder('active')        whether OPT.lsqnonlin should attach the output function
%   stop = phcpf_recorder('fn', x, ov, state)   lsqnonlin OutputFcn
%   [X, iters, resnorm] = phcpf_recorder('get')  X is a cell of single-precision iterates
persistent X iters resnorm active
if isempty(active), active = false; end
switch cmd
    case 'start'
        X = {}; iters = []; resnorm = []; active = true;
    case 'stop'
        active = false;
    case 'active'
        varargout = {active};
    case 'fn'
        x = varargin{1}; ov = varargin{2}; state = varargin{3};
        if active && any(strcmp(state, {'init', 'iter'}))
            X{end+1} = single(x);
            iters(end+1) = numel(iters);          % cumulative across warm-up solves
            if isfield(ov, 'resnorm'), resnorm(end+1) = ov.resnorm; else, resnorm(end+1) = NaN; end
        end
        varargout = {false};
    case 'get'
        varargout = {X, iters, resnorm};
    otherwise
        error('phcpf_recorder: unknown command %s', cmd);
end
end
