function poshidexp = stable_response(input)
cutoff = log(realmin('double'));
k=1;
    input(input*k>-cutoff) = -cutoff/k;
    input(input*k<cutoff) = cutoff/k;  %%% this lines is for preventing NaN values

poshidexp=input;
end