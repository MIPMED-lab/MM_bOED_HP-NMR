
using Calculus
using ScikitLearn.GridSearch
using ScikitLearnBase
using GaussianMixtures


## Entropy estimation
function H_est(comps, weis, meaas, covas)

    global mvMean = meaas;
    global mvCovar = covas;
    global mvPro = weis;
    
    # Zero-order Taylor series result for the approximation of the posterior entropy
    ZO_Expan = ZOTSE(meaas, covas, weis);

    # Second-order Taylor series result term for the approximation of the posterior entropy
    SO_Expan = SOTSE(meaas, covas, weis);

    # Aproximation of the Entropy for the posterior
    H_posterior = ZO_Expan-SO_Expan;

    return(H_posterior)

end

## Calculation of upper bound for the posterior Entropy

function H_Upper(w,E)

    comp = length(w);
    dims = size(E[1])[1];
    Hu = zeros(1,comp);

    for i in 1:comp
        hu = w[i]*(-log(w[i])+0.5*log(det(Matrix(E[i]))*(2*pi*exp(1))^dims))
        Hu[1,i] = hu;
    end

    return(sum(Hu))
end

## PDF of a multivariate Gaussian distribution

function mvGauss(x, MU, E)

    gx = (1/sqrt(det(Matrix(E))*(2*pi)^2))*exp((-1/2)*(x-MU)'*inv(Matrix(E))*(x-MU))

    return(gx)

end

## Lower bound for the posterior entropy
function H_Lower(w, E, MU)

    comp = length(w);
    Hl = zeros(1,comp);

    for i in 1:comp
        inter = zeros(1,comp);
        for j in 1:comp
            insu = w[j]*mvGauss(MU[i,:], MU[j,:], E[i]+E[j]);
            inter[1,j] = insu;
        end
        Hl[1,i] = w[i]*log(sum(inter))
    end

    return(-sum(Hl))

end

## Function for a multivariate PDF of a Gaussian Mixture
function GaussMix(x, MU, E, w)

    comp = length(w);
    FX = zeros(1,comp);

    for i in 1:comp
        fx = w[i]*(1/sqrt(det(Matrix(E[i]))*(2*pi)^2))*exp((-1/2)*(x-MU[i,:])'*inv(Matrix(E[i]))*(x-MU[i,:]));
        FX[i] = fx;
    end

    return(sum(FX))

end

## Zero-order Taylor Series Expansion

function ZOTSE(MU, E, w)

    comp = length(w);
    ZO = zeros(1,comp);

    for i in 1:comp
        zo = w[i]*log(GaussMix(MU[i,:], MU, E, w));
        ZO[i] = zo;
    end

    return(-sum(ZO))
end

## Function to calculate the most computationaly expensive part of the Taylor Expansion component
function GaussMix2(x)

    comp = length(mvPro);
    FX = zeros(1,comp);

    for i in 1:comp
        fx = mvPro[i]*(1/sqrt(det(Matrix(mvCovar[i]))*(2*pi)^2))*exp((-1/2)*(x-mvMean[i,:])'*inv(Matrix(mvCovar[i]))*(x-mvMean[i,:]));
        FX[i] = fx;
    end

    return(sum(FX))

end

function FMix(x, MU, E, w)

    comp = length(w);
    dims = size(E[1])[1];

    F_x = Array{Any,1}(undef,comp);

    for j in 1:comp
        inter = (1/GaussMix(x, MU, E, w))*(x-MU[j,:])*(Calculus.gradient(GaussMix2, x))'+(x-MU[j,:])*(inv(Matrix(E[j]))*(x-MU[j,:]))'-Diagonal(ones(dims));
        fx = w[j]*inv(Matrix(E[j]))*(inter)*(mvGauss(x , MU[j,:], E[j]));
        F_x[j] = fx;
    end
    enn = sum(F_x)/GaussMix(x, MU, E, w)

    return(enn)
end

## Calculation of the Second-Order Taylor Series expanssion term

function SOTSE(MU, E, w)

    comp = length(w);
    dims = size(E[1])[1];

    H2 = zeros(1,comp);

    for i in 1:comp
        println(string("              Iteration ", i, " of ", comp))
        inter2 = (w[i]/2)*sum(FMix(MU[i,:],MU,E,w).*Matrix(E[i]))
        H2[i] = inter2;
    end

    return(sum(H2))

end



function entroSing(matDis)
    grid_searchM1FO3 = ScikitLearnBase.fit!(GridSearchCV(GMM(n_components=10, kind=:full), 
        Dict(:n_components=>collect(1:10)), verbose=false), matDis)

    grid_searchM1FbO3 = ScikitLearnBase.fit!(GridSearchCV(GMM(n_components=grid_searchM1FO3.best_estimator_.n*4, kind=:full), 
            Dict(:n_components=>collect(grid_searchM1FO3.best_estimator_.n*4:grid_searchM1FO3.best_estimator_.n*4)), verbose=false), matDis)

    # First dimension number of experiments, second dimension type of experiments.
    comps = Array{Any,2}(undef,1, 1);
    weis = Array{Any,2}(undef, 1, 1);
    meaas = Array{Any,2}(undef,1, 1);
    covas = Array{Any,2}(undef,1, 1);

    for i in 1:1 #number of experiments
        # Get and store best results
        comps[i,1] = length(grid_searchM1FbO3.best_estimator_.w); # Number of components
        weis[i,1] = grid_searchM1FbO3.best_estimator_.w; # Weights for each component
        meaas[i,1] = grid_searchM1FbO3.best_estimator_.μ; # Means for each component
        covas[i,1] = [inv(Matrix(grid_searchM1FbO3.best_estimator_.Σ[j]))*inv(Matrix(grid_searchM1FbO3.best_estimator_.Σ[j])') 
            for j in 1:length(grid_searchM1FbO3.best_estimator_.Σ)]; # Covariance Matrix for each component
    end

    H_U_Post = Array{Any,2}(undef,1,1);
    H_L_Post = Array{Any,2}(undef,1,1);
    H_Post = Array{Any,2}(undef,1,1);

    for j in 1:1# number of experiments
        [H_U_Post[j,i] = H_Upper(weis[j,i], covas[j,i]) for i in 1:1];
        [H_L_Post[j,i] = H_Lower(weis[j,i], covas[j,i], meaas[j,i]) for i in 1:1];
        [H_Post[j,i] = H_est(comps[j,i], weis[j,i], meaas[j,i], covas[j,i]) for i in 1:1];
    end

    return(H_U_Post, H_Post, H_L_Post)
end