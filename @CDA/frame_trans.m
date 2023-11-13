function out = frame_trans(~,vec,shift,angle)                               % Linear tranformation of coordinate frames
    out = zeros(size(vec));
    for i = 1:length(vec(1,:))
        out(:,i) = [cos(angle),-sin(angle); sin(angle),cos(angle)]*vec(:,i)+shift;
    end
end