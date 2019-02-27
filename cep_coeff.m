function CC = cep_coeff(feats,M,N,L)

% M=40, N=13; L = 22;
% Type III DCT matrix routine (see Eq. (5.14) on p.77 of [1])
    dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
                                       .* repmat(pi*([1:M]-0.5)/M,N,1) ) );

    % Cepstral lifter routine (see Eq. (5.12) on p.75 of [1])
       ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
    % DCT matrix computation
    DCT = dctm( N, M );

    % Conversion of logFBEs to cepstral coefficients through DCT
    CC =  DCT * feats;

    % Cepstral lifter computation
      lifter = ceplifter( N, L );
% 
%     % Cepstral liftering gives liftered cepstral coefficients
      CC = diag( lifter ) * CC; % ~ HTK's MFCCs

end