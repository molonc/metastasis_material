function [k,x,DKL] = seedingRate(y,m)
% Implementation of the MLE seeding inference method given by Proposition 7
% and Corollary 7.4 from the supplement of the following reference:
%
% Heyde A, Reiter JG, Naxerova K, Nowak MA (2019).
% Consecutive seeding and transfer of genetic diversity in metastasis.
% PNAS 116(28):14129-37.
%
% Inputs:
%     y = matrix of clone frequencies
%         Each row is a different tumor; each column is a different clone.
%         Each row will be normalized so that it sums to 100%.
%     m = row index of primary tumor
%         Should be an integer between 1 and the number of tumors.
%         If the primary tumor is not known or included, set m=0 instead.
% Outputs:
%     k = column vector of inferred seeding rates
%         i.e. the average number of cells that migrate from the primary
%              tumor to each metastastis per cell division time
%     x = row vector of inferred source clone frequencies
%         If m>0, x is set to the primary tumor clone frequencies.
%     DKL = column vector of the KL divergence from each tumor to the
%         primary tumor (or if not designated, to the inferred source)

    % Define the number of tumors M and the number of clones N.
    [M,N] = size(y); M1 = M+1;
    
    % Normalize the clone frequency data so each tumor sums to 100%.
    y = y./repmat(sum(y,2),1,N);
    
    if ismember(m,1:M) % If the primary tumor is designated, use Proposition 7.
        x = y(m,:); % Source clone frequencies match the primary tumor.
        k(m,1) = NaN; % Seeding rate from the primary tumor to itself is not defined.
        for i = [1:(m-1),(m+1):M] % For each other tumor, calculate the seeding rate.
            fun = @(k)norm(dot(x,log(y(i,:))-(psi(abs(k)*x)-psi(abs(k)))));
            k(i,1) = abs(fminsearch(fun,1));
        end
        DKL = sum(x.*log(x./y),2); % Calculate the KL divergence.
        
    else % If the primary tumor is not designated, use Corollary 7.4.
        fun = @(kx)norm([kx(M1:end)*(log(y)-(psi(kx(1:M)'*kx(M1:end))-psi(repmat(kx(1:M)',[1 N]))))'...
                            kx(1:M)*(log(y)-(psi(kx(1:M)'*kx(M1:end))-psi(repmat(kx(1:M)',[1 N]))))...
                            sum(kx(M1:end))]);
        opt = optimset('MaxFunEvals',10^4,'MaxIter',10^4); % Set permitted number of iterations.
        kx = abs(fminsearch(fun,[ones(1,M),ones(1,N)/N],opt));
        k = kx(1:M)'; % Extract seeding rates.
        x = kx(M1:end)/sum(kx(M1:end)); % Extract source clone frequencies. 
        DKL = sum(x.*log(x./y),2); % Calculate the KL divergence.
    end    
end
