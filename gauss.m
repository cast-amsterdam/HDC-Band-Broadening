function pdf = gauss(x, h, u, s)
%Produce a peak according to the gaussian distribution.
pdf = (h / (s * sqrt(2 * pi()))) * exp(-0.5 * ((x - u) .^2) / (s .^ 2)); 
end