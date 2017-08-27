function select_vech = vec_to_vech(d)
% -------------------------------------------------
% Returns selection vector select_vech such that if
% A = vec(Psi) and B = vech(Psi)
% then B = A(select_vech)
%
% Inputs:
% d: dimension of Psi
% Outputs:
% select_vech: selection vector
% ----------------------------------------------

select_vech = find(tril(ones(d)));

end