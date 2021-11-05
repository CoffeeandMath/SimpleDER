function qout = qgen(r,v,varargin)

if max(size(varargin)) == 0
    qout = [r;v];
else
    qout = [r;v;varargin{1};varargin{2}];
end

end

