function[c, d] = integer_division(a, b)
% function[c, d] = integer_division(a, b)
% Perform integer/Euclidian division:
% a = c*b + d

c = floor(a./b);
d = mod(a, b);
