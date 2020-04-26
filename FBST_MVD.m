function evidence = FBST_MVD(X, etha, NptsInt, intFlag)
% 
% Compute evidence value against the hipotesis of the mean sample X be equal etha
%   Use:
%       evidence = FBST_MVD(X, etha, NptsInt, intFlag)
% Example:
%   X = randn(300,1) + 0.25;
%   evidence = FBST_MVD(X, 0.25);
%   evidence = FBST_MVD(X, 0.25,5000,'trap');   
%   evidence = FBST_MVD(X, 0.25,[],'3/8');
%
% Inputs:
%   X: one-dimensional sample that you want to test the average value for.
%   etha: test value of the average. If etha is a vector, test against all etha values.
%   NptsInt: Number of points used to perform the integration. Default value 500 points if not specified.
%   intFlag: indicator of the integration method to be used. Options: 'trap' for trapezoidal, '1/3' (default) for Simpsom 1/3 method, and '3/8' for Simpsom 3/8 method.
% Outputs:
%   evidence: Value (or vector) of the value of the evidence against the tested hypothesis.
% Adelino P. Silva
%
% Observation: Depends on Lambert W function. See: www.mathworks.com/help/symbolic/lambertw.html
% Is possible use the function lambertw on this repository. Thanks Cleve Moler https://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-function/
%
% Adelino P. Silva: adelinocpp@yahoo.com
%
%@Article{Silva2020,
%  Author    = {Adelino Pinheiro Silva And Maur{\'I}Lio Nunes Vieira And Adriano Vilela Barbosa},
%  Title     = {Forensic Speaker Comparison Using Evidence Interval In Full Bayesian Significance Test},
%  Journal   = {Matheamtical Problems In Enginner},
%  Year      = {2020},
%  Doi       = {10.1155/2020/2914942},
%}
% 
% Calcula a evidência contra a hipótese do valor da média da amostra X ser igual ao valor de etha
%
% Entradas:
%   X: amostra unidimensional que deseja-se testar o valor da média.
%   etha: valor de teste da média. Se etha for um vetor realiza o teste contra todos os valores de etha.
%   NptsInt: Número de pontos utilizado para pa realizar a integração. Valor padrão 500 pontos se não especificado.
%   intFlag: indicador do método de integração para ser utilizado. Opções: 'trap' para trapezoidal, '1/3' (padrão) para o método de 1/3 de Simpsom, e '3/8' para o método de 3/8 de Simpsom.
%
% Saída:
%   evidence: Valor (ou vetor) do valor da evidẽncia contra a hipótese testada.
% 
% Adelino P. Silva: adelinocpp@yahoo.com
%
% Observação: depnte da função Lambert W. Veja em: www.mathworks.com/help/symbolic/lambertw.html
%
% Por favor cite este trabalho


boolVector = 0;
if (~isscalar(etha) && isvector(etha))
    boolVector = 1;
end
if (nargin < 4) || isempty(intFlag)
    intFlag = '1/3';
end
if (nargin < 3) || isempty(NptsInt)
    NptsInt = 500;
end

% Cálculo do FBST para média de X = etha com variancia desconhecida

% =========================================================================
% === Médias com variancia desconhecida ===================================
% =========================================================================

if (boolVector == 0)
    n = length(X);
    Q = (n-1)*var(X);
    c_lin = (Q^(0.5*(n-1)))*sqrt(n)/ ( (2^(0.5*n)) *gamma(0.5)*gamma(0.5*(n-1)) );
    
    NormalGammaParam = struct;
    % -------------------------------------------------------------------------
    NormalGammaParam.Xmed    = mean(X);
    NormalGammaParam.Q       = Q;
    NormalGammaParam.s       = var(X);
    NormalGammaParam.n       = n;
    NormalGammaParam.c_lin   = c_lin;
    NormalGammaParam.RhoS    = (n-2)/Q;
    NormalGammaParam.ErrPad  = sqrt(NormalGammaParam.s/NormalGammaParam.n);
    % -------------------------------------------------------------------------
    % Calculo dos limites da função Normal_Gamma
    % -------------------------------------------------------------------------
    Xmed    = NormalGammaParam.Xmed;
    Rho_A = (n-2)/(Q*(1 + (n/Q)*(etha - Xmed)^2) );
    NormalGammaParam.Rho_A = Rho_A;
    NormalGammaParam.Rho_B = Rho_A;
    alpha   = -(n-2)/Q;
    beta    =  ((n-2)/Q)*log(Rho_A) - (Rho_A)*(1 + (n/Q)*(etha - Xmed)^2);
    exp_lambertW = exp(- beta/alpha)/alpha;
    Rho_C = exp( - (alpha*lambertw(exp_lambertW) + beta)/alpha );
    Rho_D = exp( - (alpha*lambertw(-1,exp_lambertW) + beta)/alpha );
    NormalGammaParam.Rho_C = abs(double(Rho_C));
    NormalGammaParam.Rho_D = abs(double(Rho_D));
    % -------------------------------------------------------------------------
    % mu      = etha;
    % rho     = NormalGammaParam.Rho_A;
    % log_Pn_Star = 0.5*(n-2)*log(rho) - 0.5*rho*Q*(1+ (n/Q)*(mu - mX)^2);
    % NormalGammaParam.Log_Ps_C = log_Pn_Star;
    % -------------------------------------------------------------------------
    %NptsInt = 1000;
    % -------------------------------------------------------------------------
    Rho_Ini = NormalGammaParam.Rho_C;
    Rho_Fim = NormalGammaParam.Rho_D;
    vecCRho = linspace(Rho_Ini,Rho_Fim,NptsInt);
    nu_Rho  = 0.5*(n-2)*log(vecCRho./Rho_A) -(Q/2)*(vecCRho - Rho_A) + 0.5*Rho_A*n*(etha - Xmed)^2;
    ag = 0.5*(n-1);
    bg = 2/Q;
    P_12 = gampdf(vecCRho,ag,bg);
    P3 = erf(real(sqrt(nu_Rho)));
    F = (P_12).*P3;
    I = simpsom(vecCRho,F,intFlag);
    evidence = I;
else
    evidence = zeros(length(etha),1);
    n = length(X);
    Q = (n-1)*var(X);
    c_lin = (Q^(0.5*(n-1)))*sqrt(n)/ ( (2^(0.5*n)) *gamma(0.5)*gamma(0.5*(n-1)) );
    NormalGammaParam = struct;
    % -------------------------------------------------------------------------
    NormalGammaParam.Xmed    = mean(X);
    NormalGammaParam.Q       = Q;
    NormalGammaParam.s       = var(X);
    NormalGammaParam.n       = n;
    NormalGammaParam.c_lin   = c_lin;
    NormalGammaParam.RhoS    = (n-2)/Q;
    NormalGammaParam.ErrPad  = sqrt(NormalGammaParam.s/NormalGammaParam.n);
    Xmed    = NormalGammaParam.Xmed;
    alpha   = -(n-2)/Q;
    for k = 1:length(etha)
        % -------------------------------------------------------------------------
        % Calculo dos limites da função Normal_Gamma
        % -------------------------------------------------------------------------    
        Rho_A = (n-2)/(Q*(1 + (n/Q)*(etha(k) - Xmed)^2) );
        NormalGammaParam.Rho_A = Rho_A;
        NormalGammaParam.Rho_B = Rho_A;
        
        beta    =  ((n-2)/Q)*log(Rho_A) - (Rho_A)*(1 + (n/Q)*(etha(k) - Xmed)^2);
        exp_lambertW = exp(- beta/alpha)/alpha;
        Rho_C = exp( - (alpha*lambertw(exp_lambertW) + beta)/alpha );
        Rho_D = exp( - (alpha*lambertw(-1,exp_lambertW) + beta)/alpha );
        NormalGammaParam.Rho_C = abs(double(Rho_C));
        NormalGammaParam.Rho_D = abs(double(Rho_D));
        % -------------------------------------------------------------------------
        % mu      = etha;
        % rho     = NormalGammaParam.Rho_A;
        % log_Pn_Star = 0.5*(n-2)*log(rho) - 0.5*rho*Q*(1+ (n/Q)*(mu - mX)^2);
        % NormalGammaParam.Log_Ps_C = log_Pn_Star;
        % -------------------------------------------------------------------------
        %NptsInt = 1000;
        % -------------------------------------------------------------------------
        Rho_Ini = NormalGammaParam.Rho_C;
        Rho_Fim = NormalGammaParam.Rho_D;
        vecCRho = linspace(Rho_Ini,Rho_Fim,NptsInt);
        nu_Rho  = 0.5*(n-2)*log(vecCRho./Rho_A) -(Q/2)*(vecCRho - Rho_A) + 0.5*Rho_A*n*(etha(k) - Xmed)^2;
        ag = 0.5*(n-1);
        bg = 2/Q;
        P_12 = gampdf(vecCRho,ag,bg);
        P3 = erf(real(sqrt(nu_Rho)));
        F = (P_12).*P3;
        I = simpsom(vecCRho,F,intFlag);
        evidence(k) = I;
    end,
end,
% -------------------------------------------------------------------------
