


function CompRepUtil_Entro_OEDmc(ins)

    # Input converted in the write format
    Pyr = round(ins[1], digits = 1)*1000;

    # Simulate the model for each theta and Y0 sampled
    SimulsAll = Array{Any}(undef, 1000, 3);
    for i in 1:1000
        ivss = [Pyr, 0, 0, sampsY0[i,1], sampsY0[i,2], sampsY0[i,3], 0, 0, 0, 0, 0, Pyr, 0, Pyr*(sampsTh[1,end]*0.12)];
        SimulsAll[i,1], SimulsAll[i,2], SimulsAll[i,3]  = PyruvateHP_NMR_SolveAllCp(ts, vcat(sampsTh[i,:], 0), ivss, samps);
    end

    # Define obervable matrix
    Obs = zeros(length(samps), 1000)
    [Obs[:,i]=SimulsAll[1:1000,1][i][:,13,1] for i in 1:1000];

    # Fit each time point to a distribution to acount for shape using Entropy
    EntObs = zeros(length(samps)); 

    dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    names = ["timePoint"];

    for j in 2:length(samps)
        fitts1 = Dict(); 
        bestfit1 = Dict(); 
        bestfitInd1 = zeros(1,1); 

        try
            fitts1[names[1]] = fit.(dists, Ref(Obs[j,:]));
            bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2];
            bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]];
        catch
            fitts1[names[1]] = fit.(dists2, Ref(Obs[j,:]));
            bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]+1;
            bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]];
        end

        EntObs[j] = entropy(bestfit1[names[1]]); 
    end

    HES = zeros(1,1);
    HES[1,1] = sum(EntObs)

    return(HES[1,1])
end


function TranspUtil_Entro_OEDmc(ins)

    # Input converted in the write format
    Pyr = round(ins[1], digits = 1)*1000;

    # Simulate the model for each theta and Y0 sampled
    ivss = [Pyr, 0, 0];
    SimulsAll, SimulsAll, SimulsAll  = PyruvateHP_NMR_SolveAllTb(ts, hcat(sampsTh, zeros(1000)), ivss, samps);
    
    # Define obervable matrix
    Obs = zeros(length(samps), 1000)
    [Obs[:,i]=SimulsAll[i][:,4,1] for i in 1:1000];

    # Fit each time point to a distribution to acount for shape using Entropy
    EntObs = zeros(length(samps)); 

    dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    names = ["timePoint"];

    for j in 2:length(samps)
        fitts1 = Dict(); 
        bestfit1 = Dict(); 
        bestfitInd1 = zeros(1,1); 

        try
            fitts1[names[1]] = fit.(dists, Ref(Obs[j,:]));
            bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2];
            bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]];
        catch
            fitts1[names[1]] = fit.(dists2, Ref(Obs[j,:]));
            bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]+1;
            bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]];
        end

        EntObs[j] = entropy(bestfit1[names[1]]); 
    end

    HES = zeros(1,1);
    HES[1,1] = sum(EntObs)

    return(HES[1,1])
end



function BhattacharyyaDist(mu1, mu2, sd1, sd2)
    
    E = (sd1+sd2)/2;
    Em1 = inv(E);
    dE = abs(det(E));
    
    t1 = mu1'-mu2';
    t2 = mu1-mu2;
    
    ft = (1/8)*t1*Em1*t2;
    st = dE/sqrt(abs(det(sd1))*abs(det(sd2)));
    
    bhd = ft+0.5*log(st);
    
    return(bhd);
    
end


function Utility_OEDms(ins)

    Pyr = round(ins[1], digits = 1)*1000;
    smps = size(sampsTh_CR)[1]
    # Simulate the models for each theta and Y0 sampled
    # Competitivve Repression
    SimulsAll_CR = Array{Any}(undef, smps, 3);
    for i in 1:smps
        ivss = [Pyr, 0, 0, sampsY0_CR[i,1], sampsY0_CR[i,2], sampsY0_CR[i,3], 0, 0, 0, 0, 0, Pyr, 0, Pyr*(sampsTh_CR[1,end]*0.12)];
        SimulsAll_CR[i,1], SimulsAll_CR[i,2], SimulsAll_CR[i,3]  = PyruvateHP_NMR_SolveAllCp(ts, vcat(sampsTh_CR[i,:], 0), ivss, samps);
    end

    SimulsAll_AL = Array{Any}(undef, smps, 3);
    for i in 1:smps
        ivss = [Pyr, 0, 0, sampsY0_AL[i,1], 0, 0, 0, 0, 0, Pyr, 0, Pyr*(sampsTh_AL[1,end]*0.12)];
        SimulsAll_AL[i,1], SimulsAll_AL[i,2], SimulsAll_AL[i,3]  = PyruvateHP_NMR_SolveAllAl(ts, vcat(sampsTh_AL[i,:], 0), ivss, samps);
    end

    # Define obervable matrix
    Obs_CR = zeros(length(samps)-1, smps)
    [Obs_CR[:,i]=SimulsAll_CR[1:smps,1][i][2:end,13,1] for i in 1:smps];

    Obs_AL = zeros(length(samps)-1, smps)
    [Obs_AL[:,i]=SimulsAll_AL[1:smps,1][i][2:end,11,1] for i in 1:smps];


    # Compute means and covariances
    muCR = mean!(ones(length(samps)-1), Obs_CR);
    sdCR = cov(transpose(Obs_CR)).+0.1*Matrix((Diagonal(ones(length(samps)-1))));

    muAL = mean!(ones(length(samps)-1), Obs_AL);
    sdAL = cov(transpose(Obs_AL)).+0.1*Matrix((Diagonal(ones(length(samps)-1))));

    bd = BhattacharyyaDist(muCR, muAL, sdCR, sdAL)

    return(bd)

end



function ObjectFunctMECp_MLE(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    pD2 = vcat(p[1:end-3], 0);

    # Define time vector
    t2cor = dat["sts"];

    # Define equaly-spaced time vector
    ts1 = collect(0:t2cor[end]);

    # Define initial value for simulation (use of experimental mean)
    ivss1 = [dat["Means"][1,1,1]/(ScFm*0.12), 0, dat["Means"][1,1,2]/(ScFm*0.12), p[end-2], p[end-1], p[end], 0, 0, 0, 0, 0, dat["Means"][1,1,1]/(ScFm*0.12), dat["Means"][1,1,2], dat["Means"][1,1,1]];

    # Convert sampling vector to integer to extract correct elements from simulation
    samps1 = convert.(Int, t2cor);

    # Simulate
    SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllCp(ts1, pD2, ivss1, samps1);
    
    mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat["Erros"][:,1,2].^2) .+ ((SimOnTime1[:,13] .- dat["Means"][:,1,2]).^2)./(dat["Erros"][:,1,2].^2)));

    # obj = (mm21+mm22)*(mm11+mm12);
    if p[end-2] >= p[end-1]
        mm22 = mm22*10;
    end
    obj = -(mm22);

    
    return(obj)

end
















function genStanInitDict(samps, names, chains)
    alltog = Array{Dict{String,Any},1}(undef,chains);
    for i in 1:(chains)
        global tmp = Dict{String, Any}()
        for j in 1:length(names)
            tmp[names[j]] = samps[j,i];
        end
        alltog[i] = tmp
    end

    return(alltog)
end



function convertBoundTo2(x, bo, up)

    if bo>=up
        println("Please, give correct bounds")
        return
    end
    a1, a2, b1, b2 = minimum(x), maximum(x), bo, up

    t = b1.+ ((x.-a1).*(b2-b1))./(a2-a1)
    return(t)

end

function fitPriorSampsCR(priorsamps)
    
    priorsamps = hcat(priorsamps[:,3:12], priorsamps[:,14:16]);
    thn = 13;

    priorfit = Dict();
    if size(priorsamps)[1] == thn
        priorsamps = priorsamps';
    end

    Posterior2 = zeros(size(priorsamps)); # Temporary mapped samples to assess beta distribution with the rest
    for j in 1:thn
        Posterior2[:,j] = convertBoundTo2(priorsamps[:,j], 0.1, 1)
    end

    dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform]; # For the parameters that are not between 0 and 1 (to avoid error)

    fitts = Dict(); # Where all the distribution fits will be stored
    bestfit = Dict(); # Where the best distribution fit will be stored
    bestfitInd = zeros(1,thn); # Index for the best distribution fit for each parameter

    fittsNorm = Dict();
    bestfitNorm = Dict();
    bestfitIndNorm = zeros(1,thn);

    names = ["kin", "kpl", "kbn", "kun", "kbp", "kup", "kunE", "kbnE", "ki", "kr", "NADH", "NAD", "LDH"];

    for i in 1:length(names)
        fittsNorm[names[i]] = fit.(dists, Ref(Posterior2[:,i]));
        bestfitIndNorm[i] = findmax(loglikelihood.(fittsNorm[names[i]], Ref(Posterior2[:,i])))[2];
        bestfitNorm[names[i]] = fittsNorm[names[i]][findmax(loglikelihood.(fittsNorm[names[i]], Ref(Posterior2[:,i])))[2]];
    end

    for i in 1:length(names)
        try
            fitts[names[i]] = fit.(dists, Ref(priorsamps[:,i][priorsamps[:,i] .>= 0]));
            bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2];
            bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
        catch
            fitts[names[i]] = fit.(dists2, Ref(priorsamps[:,i][priorsamps[:,i] .>= 0]));
            bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]+1;
            bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
        end
    end

    newPri = Vector{String}(undef,thn);
    parDef = Vector{String}(undef,thn);
    transPar = Vector{String}(undef,thn);

    for i in 1:thn
        parDef[i] = string("real ", names[i], "; \n");
        if i <=10
            transPar[i] = string("theta[", i+2, "] = ", names[i], "; \n");
        else
            transPar[i] = string("inits[", i-10, "] = ", names[i], "; \n");
        end

        if bestfitInd[i] == 1 # Beta
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 1) # Map samples to appropiate range
            fitSamp = fit(dists[1], conSamp) # Fit distribution to mapped samples
            newPri[i] = string(names[i], " ~ beta(", fitSamp.α, ", ", fitSamp.β, "); \n"); # Define string with pdf for prior
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",1," - (",0,")); \n"); # Define remapping to true range of values

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",1," - (",0,")); \n");;
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",1," - (",0,")); \n");;
            end

        elseif bestfitInd[i] == 2 # Exponential
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 2)
            fitSamp = fit(dists[2], conSamp)
            newPri[i] = string(names[i], " ~ exponential(", fitSamp.θ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");
            end

        elseif bestfitInd[i] == 3 # LogNormal
            parDef[i] = string("real ", names[i], "; \n");
            # transPar[i] = string("theta[",i,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n"); # Map logNormal to univariate Normal. This is how to transform the samples back
            newPri[i] = string(names[i], " ~ normal(", 0, ", ", 1, "); \n"); # Univariate normal

            if i <=10
                transPar[i] = string("theta[",i+2,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            end

        elseif bestfitInd[i] == 4 # Normal
            if ((0-bestfit[names[i]].μ)/bestfit[names[i]].σ) >=-2 # Truncation only if the actual 0 for the parameter is inside the defined distribution.
                parDef[i] = string("real<lower=",(0-bestfit[names[i]].μ)/bestfit[names[i]].σ,"> ", names[i], "; \n");
            else
                parDef[i] = string("real ", names[i], "; \n");
            end
            # transPar[i] = string("theta[",i,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n") # Map Normal to univariate Normal. This is how to transform the samples back
            newPri[i] = string(names[i], " ~ normal(", 0, ", ", 1, "); \n"); # Univariate normal

            if i <=10
                transPar[i] = string("theta[",i+2,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            end

        elseif bestfitInd[i] == 5 # Gamma
            conSamp = convertBoundTo2(priorsamps[:,i], 0.001, 2)
            fitSamp = fit(dists[5], conSamp)
            newPri[i] = string(names[i], " ~ gamma(", fitSamp.α, ", ", fitSamp.θ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");
            end

        elseif bestfitInd[i] == 6 # Laplace
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[6], conSamp)
            newPri[i] = string(names[i], " ~ double_exponential(", fitSamp.μ, ", ", fitSamp.θ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            end

        elseif bestfitInd[i] == 7 # Pareto
            conSamp = convertBoundTo2(priorsamps[:,i], 0.001, 2)
            fitSamp = fit(dists[7], conSamp)
            newPri[i] = string(names[i], " ~ pareto(", fitSamp.θ, ", ", fitSamp.α, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");
            end

        elseif bestfitInd[i] == 8 # Rayleigh
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 2)
            fitSamp = fit(dists[8], conSamp)
            newPri[i] = string(names[i], " ~ rayleigh(", fitSamp.σ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");
            end

        elseif bestfitInd[i] == 9 # Cauchy
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[9], conSamp)
            newPri[i] = string(names[i], " ~ cauchy(", fitSamp.μ, ", ", fitSamp.σ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            end

        elseif bestfitInd[i] == 10
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[10], conSamp)
            newPri[i] = string(names[i], " ~ uniform(", -2, ", ", 2, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            end

        end
    end

    
    pushfirst!(transPar, "theta[2] = 55; \n");
    pushfirst!(transPar, "theta[1] = 48; \n");
    pushfirst!(transPar, "array[3] real inits; \n");
    pushfirst!(transPar, "real theta[nParms]; \n");
    push!(transPar, "theta[13] = 15; \n");
    
    

    priorfit["pars"] = parDef;
    priorfit["transpars"] = transPar;
    priorfit["pridis"] = newPri;

    return(priorfit)
end


function fitPriorSampsCR_SIMP(priorsamps)
    
    priorsamps = hcat(priorsamps[:,3:12], priorsamps[:,14:16]);
    thn = 13;

    priorfit = Dict();
    if size(priorsamps)[1] == thn
        priorsamps = priorsamps';
    end

    Posterior2 = zeros(size(priorsamps)); # Temporary mapped samples to assess beta distribution with the rest
    for j in 1:thn
        Posterior2[:,j] = convertBoundTo2(priorsamps[:,j], 0.1, 1)
    end

    dists = [Exponential, LogNormal, Normal, Laplace, Cauchy, Uniform];
    dists2 = [Exponential, LogNormal, Normal, Laplace, Cauchy, Uniform]; # For the parameters that are not between 0 and 1 (to avoid error)

    fitts = Dict(); # Where all the distribution fits will be stored
    bestfit = Dict(); # Where the best distribution fit will be stored
    bestfitInd = zeros(1,thn); # Index for the best distribution fit for each parameter

    fittsNorm = Dict();
    bestfitNorm = Dict();
    bestfitIndNorm = zeros(1,thn);

    names = ["kin", "kpl", "kbn", "kun", "kbp", "kup", "kunE", "kbnE", "ki", "kr", "NADH", "NAD", "LDH"];

    for i in 1:length(names)
        fittsNorm[names[i]] = fit.(dists, Ref(Posterior2[:,i]));
        bestfitIndNorm[i] = findmax(loglikelihood.(fittsNorm[names[i]], Ref(Posterior2[:,i])))[2];
        bestfitNorm[names[i]] = fittsNorm[names[i]][findmax(loglikelihood.(fittsNorm[names[i]], Ref(Posterior2[:,i])))[2]];
    end

    for i in 1:length(names)
        try
            fitts[names[i]] = fit.(dists, Ref(priorsamps[:,i][priorsamps[:,i] .>= 0]));
            bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2];
            bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
        catch
            fitts[names[i]] = fit.(dists2, Ref(priorsamps[:,i][priorsamps[:,i] .>= 0]));
            bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]+1;
            bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
        end
    end

    newPri = Vector{String}(undef,thn);
    parDef = Vector{String}(undef,thn);
    transPar = Vector{String}(undef,thn);

    for i in 1:thn
        parDef[i] = string("real ", names[i], "; \n");
        if i <=10
            transPar[i] = string("theta[", i+2, "] = ", names[i], "; \n");
        else
            transPar[i] = string("inits[", i-10, "] = ", names[i], "; \n");
        end

        if bestfitInd[i] == 1 # Exponential
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 2)
            fitSamp = fit(dists[1], conSamp)
            newPri[i] = string(names[i], " ~ exponential(", fitSamp.θ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");
            end

        elseif bestfitInd[i] == 2 # LogNormal
            parDef[i] = string("real ", names[i], "; \n");
            # transPar[i] = string("theta[",i,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n"); # Map logNormal to univariate Normal. This is how to transform the samples back
            newPri[i] = string(names[i], " ~ normal(", 0, ", ", 1, "); \n"); # Univariate normal

            if i <=10
                transPar[i] = string("theta[",i+2,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            end

        elseif bestfitInd[i] == 3 # Normal
            if ((0-bestfit[names[i]].μ)/bestfit[names[i]].σ) >=-2 # Truncation only if the actual 0 for the parameter is inside the defined distribution.
                parDef[i] = string("real<lower=",(0-bestfit[names[i]].μ)/bestfit[names[i]].σ,"> ", names[i], "; \n");
            else
                parDef[i] = string("real ", names[i], "; \n");
            end
            # transPar[i] = string("theta[",i,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n") # Map Normal to univariate Normal. This is how to transform the samples back
            newPri[i] = string(names[i], " ~ normal(", 0, ", ", 1, "); \n"); # Univariate normal

            if i <=10
                transPar[i] = string("theta[",i+2,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n");
            end

        elseif bestfitInd[i] == 4 # Laplace
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[4], conSamp)
            newPri[i] = string(names[i], " ~ double_exponential(", fitSamp.μ, ", ", fitSamp.θ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            end

        elseif bestfitInd[i] == 5 # Cauchy
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[5], conSamp)
            newPri[i] = string(names[i], " ~ cauchy(", fitSamp.μ, ", ", fitSamp.σ, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            end

        elseif bestfitInd[i] == 6
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[6], conSamp)
            newPri[i] = string(names[i], " ~ uniform(", -2, ", ", 2, "); \n");
            # transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

            if i <=10
                transPar[i] = string("theta[",i+2,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            else
                transPar[i] = string("inits[",i-10,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");
            end

        end
    end

    
    pushfirst!(transPar, "theta[2] = 55; \n");
    pushfirst!(transPar, "theta[1] = 48; \n");
    pushfirst!(transPar, "array[3] real inits; \n");
    pushfirst!(transPar, "real theta[nParms]; \n");
    push!(transPar, "theta[13] = 15; \n");
    
    

    priorfit["pars"] = parDef;
    priorfit["transpars"] = transPar;
    priorfit["pridis"] = newPri;

    return(priorfit)
end




function GenStanModel_CR(pris, iter, type)

    ModelSeq1 = string("
    
    functions{
    
       real[] PyruvateHP_NMR_ODEs4(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    
          // Parameters
          real T1_X = p[1]; 
          real T1_P = p[2];
          real kin = p[3];
          real kpl = p[4];
          real kbn = p[5];
          real kun = p[6];
          real kbp = p[7];
          real kup = p[8];
          real kunE = p[9];
          real kbnE = p[10];
          real ki = p[11];
          real kr = p[12];
          real ScFm = p[13];
    
          // ODEs
          real dInd_dt[14];
    
          real Pout = y[1];
          real Php = y[2]; 
          real Xhp = y[3]; 
          real NADH = y[4];
          real NAD = y[5];
          real LDH = y[6];
          real LDHac = y[7];
          real LDHacE = y[8];
          real LDHna = y[9];
          real LDHP = y[10];
          real XhpHyper = y[11];
          real PhpHyper = y[12];
          real Obs_XhpHyper = y[13];
          real Obs_PhpHyper = y[14];
    
          dInd_dt[1] = - (Pout*kin);
          dInd_dt[2] = (Pout*kin) + kup*LDHP + kr*LDHna - kbp*LDHac*Php - ki*LDH*Php; 
          dInd_dt[3] = kpl*LDHP; 
          dInd_dt[4] = kun*LDHac - kbn*LDH*NADH;
          dInd_dt[5] = kunE*LDHacE - kbnE*LDH*NAD;
          dInd_dt[6] = kun*LDHac + kr*LDHna + kunE*LDHacE - kbn*LDH*NADH - ki*LDH*Php - kbnE*LDH*NAD;
          dInd_dt[7] = kbn*LDH*NADH + kup*LDHP - kun*LDHac - kbp*LDHac*Php;
          dInd_dt[8] = kpl*LDHP + kbnE*LDH*NAD - kunE*LDHacE;
          dInd_dt[9] = ki*LDH*Php - kr*LDHna;
          dInd_dt[10] = kbp*LDHac*Php - kup*LDHP - kpl*LDHP;
          dInd_dt[11] = dInd_dt[3] - (XhpHyper/T1_X);
          dInd_dt[12] = (dInd_dt[1] + dInd_dt[2] + dInd_dt[9] + dInd_dt[10]) - (PhpHyper/T1_P);
          dInd_dt[13] = dInd_dt[11] * (ScFm*0.12);
          dInd_dt[14] = dInd_dt[12] * (ScFm*0.12);
    
          // Results
          return dInd_dt;
    
        }
    
        vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    
          return init; 
    
        }
    
    }
    
    
    data {
        // Observables
        int m; // Total number of data series
        int stslm; // Maximum number of rows for all the observable matrices
        array[1,m] int stsl; // Number of elements at each time series for each series m
        array[stslm,m] int sts; // Sampling times for each series m
    
        int obser;//-> Introduce this so we have all the data in one same array (easier generalisation). Work on generalisation in case different experiments have different obsevables?
        array[1,obser] int obSta; // -> This variable will be to know which are the observable states
        array[1,m] int ncells;
        array[1,m] int nts;
    
       //  real itp[m];
        array[stslm,m,obser] real Means; // ---> General arrays of means and errors
        array[stslm,m,obser] real Erros;
    
        // Inputs
        int tml; // Maximum length of the rows for the sampling times matrix
        array[tml, m] real ts; // Time series for each serie m
        array[14,m] real Y0us; // Y0 vectors
    }
    
    
    transformed data {
        int nParms = 13; // Number of parameters of the model //--------> Introduce number in generation of script
        int Neq = 14; // Total number of equations of the model //-----> Introduce number in generation of script
        array[0] int x_i; // Empty x_i object (needs to be defined)
        array[0] real x_r; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
        array[Neq,m] real ivss = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
        real itp = 0;   
    }
    
    
    parameters {
    
        ",join(pris["pars"]),"
    }
    
    
    transformed parameters {
    
       ",join(pris["transpars"]),"
    
    }
    
    
    
    
    model {
       // Intermediate parameters
      vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
      
      array[Neq,m] real Y0; // Initial values for the ODEs variables at the first event
      array[stslm,m,Neq] real yhat; // ---> Generall array to include all the observables (easier generalisation)
    
    ",join(pris["pridis"]),"
    
    
       // Likelihood
       for (j in 1:m){
    
          array[Neq] real ivst; // Initial value of the states
          array[tml+1,Neq] real y_hat; // Object to include the ODEs solutions for each state
          int lts = num_elements(ts[1:nts[1,j]+1,j]);  // Length of the time series for each event
          array[lts,Neq] real part1; // Temporary object that will include the solution of the ODE for each event at each loop
          vector[stsl[1,j]] yhatPyr;
          // Loop (over the number of events/inputs) to solve the ODE for each event stopping the solver and add them to the final object y_hat
    
          int q=1;
          
    
          // Calculation of the solution for the ODEs where for events that are not the first one. The time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
          ivst = ivss[,j];
          
          ivst[4] = inits[1];
          ivst[5] = inits[2];
          ivst[6] = inits[3];
    
          // print(itp);
          
          part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs4,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-9, 1.0e-9, 1e7);
          
          // Introduction of the result of part1 into the object y_hat
          y_hat[(1),] = ivst;
          for (y in 1:lts){
             y_hat[(y)+1,]=part1[(y),];
          };
    
          for (t in 1:stsl[1,j]){ //----> General form
             for (ob in 1:Neq){
                yhat[t,j,ob] = part1[(sts[t,j]+1),ob];
             }
          }
    
    
          
    
          // Means[1:stsl[1,j],j,1] ~ normal(yhatPyr[1:stsl[1,j]],Erros[1:stsl[1,j],j,1]);
          // target += normal_lpdf(Means[1:stsl[1,j],j,1] | yhat[1:stsl[1,j],j,14], Erros[1:stsl[1,j],j,1]);
    
    
          // Means[1:stsl[1,j],j,2] ~ normal(yhat[1:stsl[1,j],j,3],Erros[1:stsl[1,j],j,2]);
          target += normal_lpdf(Means[1:stsl[1,j],j,2] | yhat[1:stsl[1,j],j,13], Erros[1:stsl[1,j],j,2]);
    
    
       }
    }
    
    
        
        ");
    
        open(string("C:\\IBECPostDocDrive\\2024_01_16_NCvsKR\\DataProcessingInference\\ComputationalBayesianOED\\StanModels\\HP_PyrLac_CompetitiveRepressionModel_FixedT1ScF_",type,"_Iter",iter,".stan"), "w") do io
            write(io, ModelSeq1);
         end;
    
             end




function dPSD_M(poster, trueP, upB, lwB)

    dPSD = zeros(size(poster));

    for i in 1:size(poster)[1]
        dPSD[i,:] = (poster[i,:] - trueP)./(upB - lwB);
    end

    dPSDM = sum(sqrt.(dPSD.^2), dims = 1)./size(dPSD)[1];

    return dPSD, dPSDM

end

























