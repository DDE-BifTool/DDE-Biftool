function x=dde_complex_from_point(point,free_par_ind)
[ind,len]=dde_ind_from_point(point,free_par_ind);
fnames=fieldnames(ind);
point.parameter=point.parameter(free_par_ind);
z=zeros(len,1);
rs=@(x)reshape(x,[],1);
for k=1:length(fnames)
    fd=ind.(fnames{k});
    if ~isstruct(fd)
        z(rs(fd))=rs(point.(fnames{k}));
    else
        z(rs(fd.re))=rs(point.(fnames{k}));
        z(rs(fd.im))=rs(1i*point.(fnames{k}));
    end
end
x=[real(z),imag(z)];
end
