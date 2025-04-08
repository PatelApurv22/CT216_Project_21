baseGraph5GNR = 'NR_2_6_52'; % load 5G NR LDPC base H matrix, use both NR_2_6_52 and NR_1_5_352

codeRate = 1/4; % vary this in the set {1/4, 1/3, 1/2, 3/5} for 2_6_52 and in the range {1/3, 1/2, 3/5 and 4/5} for 1_5_352

[B,Hfull,z] = nrldpc_Hmatrix(baseGraph5GNR); % Convert the base H matrix to binary H matrix

[mb,nb] = size(B); kb = nb - mb; % 5G NR specific details

kNumInfoBits = kb * z; % Number of information bits

k_pc = kb-2; nbRM = ceil(k_pc/codeRate)+2; % Some 5G NR specific details 

nBlockLength = nbRM * z; % Number of encoded bits

% Next three lines are some 5G NR specific details
H = Hfull(:,1:nBlockLength); 

nChecksNotPunctured = mb*z - nb*z + nBlockLength;
H = H(1:nChecksNotPunctured,:); % this is the binary H matrix

Nchecks = size(H,1); % Number of CNs (we have denoted this as U = N - K in the class)
b = randi([0 1],[kNumInfoBits 1]); % Generate information (or message) bit vector
c = nrldpc_encode(B,z,b'); % Encode using 5G NR LDPC base matrix 
c = c(1:nBlockLength)';

function out = mul_sh(vec, sh)
    if sh == -1
        out = zeros(1, length(vec));
    else
        out = circshift(vec, [0 -sh]);
    end
end

function [B,H,z] = nrldpc_Hmatrix(BG)
    load(sprintf('%s.txt',BG),BG);
    B = NR_2_6_52;
    [mb,nb] = size(B);
    z = 52;
    H = zeros(mb*z,nb*z); 
    Iz = eye(z); I0 = zeros(z);
    for kk = 1:mb
        tmpvecR = (kk-1)*z+(1:z);
        for kk1 = 1:nb
            tmpvecC = (kk1-1)*z+(1:z);
            if B(kk,kk1) == -1
                H(tmpvecR,tmpvecC) = I0;
            else
                H(tmpvecR,tmpvecC) = circshift(Iz,-B(kk,kk1));
            end
        end
    end
    
    [U,N]=size(H); K = N-U;
    P = H(:,1:K);
    G = [eye(K); P];
    Z = H*G;

end

function cword = nrldpc_encode(B,z,msg)
    %B: base matrix
    %z: expansion factor
    %msg: message vector, length = (#cols(B)-#rows(B))*z
    %cword: codeword vector, length = #cols(B)*z
    
    [m,n] = size(B);
    
    cword = zeros(1,n*z);
    cword(1:(n-m)*z) = msg;
    
    %double-diagonal encoding
    temp = zeros(1,z);
    for i = 1:4 %row 1 to 4
        for j = 1:n-m %message columns
            temp = mod(temp + mul_sh(msg((j-1)*z+1:j*z),B(i,j)),2);
        end
    end
    if B(2,n-m+1) == -1
        p1_sh = B(3,n-m+1);
    else
        p1_sh = B(2,n-m+1);
    end
    cword((n-m)*z+1:(n-m+1)*z) = mul_sh(temp,z-p1_sh); %p1
    %Find p2, p3, p4
    for i = 1:3
        temp = zeros(1,z);
        for j = 1:n-m+i
            temp = mod(temp + mul_sh(cword((j-1)*z+1:j*z),B(i,j)),2);
        end
        cword((n-m+i)*z+1:(n-m+i+1)*z) = temp;
    end
    %Remaining parities
    for i = 5:m
        temp = zeros(1,z);
        for j = 1:n-m+4
            temp = mod(temp + mul_sh(cword((j-1)*z+1:j*z),B(i,j)),2);        
        end
        cword((n-m+i-1)*z+1:(n-m+i)*z) = temp;    
    end

end