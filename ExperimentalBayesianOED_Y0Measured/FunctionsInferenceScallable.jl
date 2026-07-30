function PyruvateHP_CellsCp_Scall!(du, u, p, t)

    Pout, Php, Xhp, NADH, NAD, LDH, LDHac, LDHacE, LDHna, LDHP, XhpHyper, PhpHyper, Obs_XhpHyper, Obs_PhpHyper = u;
    T1_X, T1_P, kin, kpl, kbn, kun, kbp, kup, kunE, kbnE, ki, kr, ClNm, ScFm = p;


    # kin = kin*ClNm
    # Xhp, NADH, NAD, LDH, LDHac, LDHacE, LDHna, LDHP = Xhp*ClNm, NADH*ClNm, NAD*ClNm, LDH*ClNm, LDHac*ClNm, LDHacE*ClNm, LDHna*ClNm, LDHP*ClNm



    du[1] = dPout = - (Pout*kin)
    du[2] = dPhp = (Pout*kin) + kup*LDHP + kr*LDHna - kbp*LDHac*Php - ki*LDH*Php
    du[3] = dXhp = kpl*LDHP
    du[4] = dNADH = kun*LDHac - kbn*LDH*NADH
    du[5] = dNAD = kunE*LDHacE - kbnE*LDH*NAD
    du[6] = dLDH = kun*LDHac + kr*LDHna + kunE*LDHacE - kbn*LDH*NADH - ki*LDH*Php - kbnE*LDH*NAD
    du[7] = dLDHac = kbn*LDH*NADH + kup*LDHP - kun*LDHac - kbp*LDHac*Php
    du[8] = dLDHacE = kpl*LDHP + kbnE*LDH*NAD - kunE*LDHacE
    du[9] = dLDHna = ki*LDH*Php - kr*LDHna
    du[10] = dLDHP = kbp*LDHac*Php - kup*LDHP - kpl*LDHP
    

    # du[9] = dPHyper =  dPout + dPhp - (PHyper/T1_P)
    du[11] = dXhpHyper =  dXhp*ClNm - (XhpHyper/T1_X)
    du[12] = dPhpHyper =  dPout + dPhp + dLDHna + dLDHP - (PhpHyper/T1_P)

    du[13] = dObs_XhpHyper =  dXhpHyper * (ScFm*0.12)
    du[14] = dObs_PhpHyper =  dPhpHyper * (ScFm*0.12)

end


function PyruvateHP_NMR_SolveAllCp_Scall(ts, pD, ivss, samps)


    nStat = 14;
    nPar = 15;

    if length(size(pD)) == 1
        pD = reshape(pD,size(pD)[1],1);
    end

    if size(pD)[2] != nPar     
        pD = pD';
    end

    if length(ivss)/nStat > 1
        if size(ivss)[2] != nStat
            ivss = ivss';
        end
    end

    AllSolTest = zeros(length(samps), nStat, length(pD[:,1])); # Simulation of the system observed
    AllSolTest_Off = Array{Any,1}(undef,length(pD[:,1])); # Simulation of the system before we obvserve it (considering time offset). First column is the time vector
    AllSolTest_Tog = Array{Any,1}(undef,length(pD[:,1])); # Previous two together. First column is the time vector. 


    
    for drawInd in collect(1:length(pD[:,1]))
        
        p = pD[drawInd,1:end-1];
        tau = pD[drawInd,end];
        
        if length(ivss)/nStat > 1
            ivss2 = ivss[drawInd,:];
        else
            ivss2 = ivss;
        end
    
        # Account for the time offset
    
        
            prob_off = ODEProblem(PyruvateHP_CellsCp_Scall!,ivss2,(-tau, 0),p);
            part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9);
        
            ivss2 = part1_off.u[end];
        
            prob = ODEProblem(PyruvateHP_CellsCp_Scall!,ivss2,(ts[1], ts[end]),p);
            part1 = DifferentialEquations.solve(prob, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
        
        
        AllSolTest[:,:,drawInd]
    
        tmp = zeros(length(part1_off.u), nStat);
    
        
        

        for j in 1:nStat
            AllSolTest[:,j,drawInd] = [part1.u[i][j] for i in 1:length(part1.u)][samps.+1];
            tmp[:,j] = [part1_off.u[i][j] for i in 1:length(part1_off.u)];
        end

        AllSolTest_Off[drawInd] = hcat(part1_off.t, tmp);
        AllSolTest_Tog[drawInd] = vcat(AllSolTest_Off[drawInd][1:end-1, :], hcat(samps, AllSolTest[:,:,drawInd]));
        
        
        
    
    end

    return AllSolTest, AllSolTest_Off, AllSolTest_Tog

end


###########################################################################################




function genPriorSamps(NSamps)

    priorTheta = Array{Any}(undef, 14);
    priotY0 = Array{Any}(undef, 3);

    priorTheta[1] = Normal(40.92, 10); 
    priorTheta[2] = Normal(51.23, 10); 
    priorTheta[3] = Truncated(Normal(0.0025, 0.0025), 0,2); 
    priorTheta[4] = Truncated(Normal(0.0275, 0.0275), 0,2); 
    priorTheta[5] = Truncated(Normal(5e-7, 5e-7), 0,2); 
    priorTheta[6] = Truncated(Normal(0.15, 0.15), 0,0.3); 
    priorTheta[7] = Truncated(Normal(0.95, 0.95), 0.05,5); 
    priorTheta[8] = Truncated(Normal(0.95, 0.95), 0,2); 
    priorTheta[9] = Truncated(Normal(0.1, 0.1), 0,2); 
    priorTheta[10] = Truncated(Normal(3e-6, 3e-6), 0,1.0e-5);
    priorTheta[11] = Truncated(Normal(0.95, 0.95), 0,5); 
    priorTheta[12] = Truncated(Normal(0.0009, 0.0009), 0,2); 
    
    priorTheta[13] = Normal(18.5, 2.29); 
    priorTheta[14] = Uniform(3e6-1,3e6)

    priotY0[1] = Truncated(Normal(457, 238), 45.7, 1409); 
    priotY0[2] = Truncated(Normal(3146, 1270), 31.46, 8226.67); 
    priotY0[3] = Truncated(Normal(0.3947, 0.193), 0, 1.168);


    global sampsTh = zeros(NSamps, 14);
    global sampsY0 = zeros(NSamps, 3);

    Random.seed!(687865443); 

    for i in 1:14
        sampsTh[:,i] = rand(priorTheta[i], NSamps);
    end
    for i in 1:3
        sampsY0[:,i] = rand(priotY0[i], NSamps);
    end

    sampsTh[:,2] .= 51.23; # T1 of Pyruvate is fixed since it is irrelevant for the model for the problem and can be estimated with experiments without cells. 
    sampsTh[:,13] .= 3e6; # In this part of the work, we will always use 3 million cells for the experiments
    sampsTh[:,14] .= 20; # This will be re-writen with the empirical value once extracted from data

    return priorTheta, priotY0, sampsTh, sampsY0

end





################################################################################################




function genThetaBounds()
    bounds = [zeros(13), zeros(13)];
    boundY0 = [[45.7, 31.46, 0], [1409, 8226.67, 1.168]];
    for i in 1:13
        if i >= 3 && i <13
            if priorTheta[i].untruncated.μ - (priorTheta[i].untruncated.σ*4) >= 0
                bounds[1][i] = priorTheta[i].untruncated.μ - (priorTheta[i].untruncated.σ*4);
            else
                bounds[1][i] = 0;
            end
            if priorTheta[i].untruncated.μ + (priorTheta[i].untruncated.σ*4) <= 2
                bounds[2][i] = priorTheta[i].untruncated.μ + (priorTheta[i].untruncated.σ*4);
            else
                bounds[2][i] = 2;
            end
        elseif i == 1
            bounds[1][i] = priorTheta[i].μ - (priorTheta[i].σ*2);
            bounds[2][i] = priorTheta[i].μ + (priorTheta[i].σ*4);
        elseif i == 2
            bounds[1][i] = 51.23;
            bounds[2][i] = 51.23;
        end
        if i == 13
            bounds[1][i] = 3e6;
            bounds[2][i] = 3e6;
        end
        if i == 10
            bounds[2][i] = 1.0e-5;
        end
        if i == 6
            bounds[2][i] = 0.3;
        end
        if i == 7
            bounds[1][i] = 0.05;
        end
    end

    return bounds, boundY0

end




###############################################################################################



function restructDatInfCompRep_RealExpDat_Scal(paths, conce, NADHmean, NADmean, LDHmean, clNm)

    m = length(paths);

    DataAll = Array{Any,1}(undef, m)

    for i in 1:m
        DataAll[i] = Matrix(CSV.read(paths[i], DataFrame));
    end
    stslm = maximum([size(DataAll[i])[1] for i in 1:m]);
    stsl = reshape([size(DataAll[i])[1] for i in 1:m], 1, m);

    obser = 2;
    obSta = Array{Int,2}(undef, 1, 2)
    obSta[1,:] = [1,2];

    Means = Array{Float64,3}(undef, stslm, m, obser);
    Erros = Array{Float64,3}(undef, stslm, m, obser);
    sts = Array{Int,2}(undef, stslm, m).*0
    sts2 = Array{Int,2}(undef, stslm+1, m)

    Y0us = Array{Float64,2}(undef, 14, m);
    tml = maximum([length(DataAll[i][1,1]:1:DataAll[i][end,1]) for i in 1:m]);
    ts = Array{Float64,2}(undef, tml, m);
    itp = Array{Float64,1}(undef, m);

    ms = collect(2:2:obser*2);

    ScFms = zeros(1, m);


    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        ScFms[i] = DataAll[i][1,2]./((conce[i]*1000)*0.12)[1];
        Y0us[:,i] = [DataAll[i][1,2]/(ScFms[i]*0.12), 0,0, NADHmean, NADmean, LDHmean, 0, 0, 0, 0, 0, DataAll[i][1,2]/(ScFms[i]*0.12), 0, DataAll[i][1,2]]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
            ncells[1,i] = 3;
    end
    ncells = convert.(Int, ncells);

    nts = zeros(1,m);
    for i in 1:m
        nts[1,i] = round(DataAll[i][end,1])
    end
    nts = convert.(Int, nts);

    # if length(ncells) == 1
    #     ncells = ncells[1];
    # end

    nmCls = zeros(1, m);
    if length(clNm) == 1
        nmCls .= clNm
    else
        try
            for i in 1:length(clNm)
                nmCls[i] = clNm
            end
        catch
            println("If different number of cells in each experiment, indicate the cell number for each experiment loaded. You did not, so this will generate and Error!")
        end
    end


    data_multi = Dict(
        
        "m" => m,
        "stslm" => stslm,
        "stsl" => stsl,
        "sts" => sts,
        # "sts2" => sts2,
        "obser" => obser,
        "obSta" => obSta,
        "Means" => Means,
        "Erros" => Erros,
        "tml" => tml,
        "ts" => ts,
        "Y0us" => Y0us,
        "itp" => itp,
        "ncells" => ncells,
        "nts" => nts,
        "ScFms" => ScFms,
        "nmCls" => nmCls
    
    );

    return data_multi

end




#########################################################################################################


function simPriorPredictAll(NSamps, dat, sampsTh, sampsY0, conce)

    simsPrior = Array{Any}(undef, NSamps);


    for mm in 1:length(dat["ScFms"])
        sampsTh[:,14] .= dat["ScFms"][mm];

        tsC2 = dat["sts"][:,mm];

        for i in 1:NSamps

            ts = 0:tsC2[end];
            ivss = [round(conce[mm], digits = 1)*1000, 0, 0, sampsY0[i,1], sampsY0[i,2], sampsY0[i,3], 0, 0, 0, 0, 0, round(conce[mm], digits = 1)*1000, 0, round(conce[mm], digits = 1)*1000*(dat["ScFms"][mm]*0.12)];;
            samps = convert.(Int, tsC2);

            SimulsAll = []
            try
                SimulsAll, SimulsAll2, SimulsAll3 = PyruvateHP_NMR_SolveAllCp_Scall(ts, vcat(sampsTh[i,1:end], 0), ivss, samps);
            catch
                tsC2 = dat["sts"][1:end-1,mm];
                ts = 0:tsC2[end];
                samps = convert.(Int, tsC2);

                SimulsAll, SimulsAll2, SimulsAll3 = PyruvateHP_NMR_SolveAllCp_Scall(ts, vcat(sampsTh[i,1:end], 0), ivss, samps);
            end
            

            simsPrior[i] = SimulsAll;
        end


        # Dims = Up/Down, Observable, Experiment
        PriorQuant = Array{Any}(undef,2,1,1);
        lacPirorSims = zeros(size(simsPrior[1])[1], NSamps)
        [lacPirorSims[:,i] = simsPrior[i][:,13,1] for i in 1:NSamps ];

        for k in 1:1
            for m in 1:1
                PriorQuant[1,k,m] = [percentile(lacPirorSims[j,:], 99.5) for j in 1:size(lacPirorSims)[1]]; # Up
                PriorQuant[2,k,m] = [percentile(lacPirorSims[j,:], 0.5 ) for j in 1:size(lacPirorSims)[1]]; # Down
            end
        end

        pr = plot(tsC2./60, PriorQuant[1], label = "", color = "blue", linewidth = 3,
                margin=10Plots.mm,xtickfont=font(16), ytickfont=font(16), guidefont=font(18), titlefont=font(18),
                xlabel = "time (min)", ylabel = "Lac (AU)", title = string(conce[mm]))
            plot!(tsC2./60, PriorQuant[2], label = "", color = "blue", linewidth = 3)
            plot!(tsC2./60, PriorQuant[1], fillrange=PriorQuant[2], label="", color="blue3", fillalpha=0.2)

            plot!(dat["sts"]./60, dat["Means"][:,mm,2], yerror = dat["Erros"][:,mm,2], linewidth = 2, colour = "black", label = "")


        PriorQuant2 = Array{Any}(undef,2,1,1);
        pyrPirorSims = zeros(size(simsPrior[1])[1], NSamps)
        [pyrPirorSims[:,i] = simsPrior[i][:,14,1] for i in 1:NSamps ];

        for k in 1:1
            for m in 1:1
                PriorQuant2[1,k,m] = [percentile(pyrPirorSims[j,:], 99.5) for j in 1:size(pyrPirorSims)[1]]; # Up
                PriorQuant2[2,k,m] = [percentile(pyrPirorSims[j,:], 0.5 ) for j in 1:size(pyrPirorSims)[1]]; # Down
            end
        end

        pr2 = plot(tsC2[1:end]./60, PriorQuant2[1][1:end], label = "", color = "blue", linewidth = 3,
            margin=10Plots.mm,xtickfont=font(16), ytickfont=font(16), guidefont=font(18), titlefont=font(18),
            xlabel = "time (min)", ylabel = "Pyr (AU)")
          plot!(tsC2[1:end]./60, PriorQuant2[2][1:end], label = "", color = "blue", linewidth = 3)
          plot!(tsC2[1:end]./60, PriorQuant2[1][1:end], fillrange=PriorQuant2[2][1:end], label="", color="blue3", fillalpha=0.2)
        
          plot!(dat["sts"][1:end,mm]./60, dat["Means"][:,mm,1][1:end], yerror = dat["Erros"][:,mm,1][1:end], linewidth = 2, colour = "black", label = "")

        # plss[mm] = plot(pr, pr, size=(1400,500))
        display(plot(pr, pr2, size=(1400,500)))
    end
    # plot(plss...)
end




#########################################################################################################



function ObjectFunctMECp_Multi_MLE2(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    
    obj = 0;
    
    for i in 1:length(dat["stsl"])
        pD2 = vcat(p[1:end-3], dat["ScFms"][i], 0);

        # Define time vector
        t2cor = dat["sts"][:,i];

        # Define equaly-spaced time vector
        ts1 = collect(0:t2cor[end]);

        # Define initial value for simulation (use of experimental mean)
        ivss1 = [dat["Means"][1,i,1]/(dat["ScFms"][i]*0.12), 0, dat["Means"][1,i,2]/(dat["ScFms"][i]*0.12), p[end-2], p[end-1], p[end], 0, 0, 0, 0, 0, dat["Means"][1,i,1]/(dat["ScFms"][i]*0.12), dat["Means"][1,i,2], dat["Means"][1,i,1]];

        
        # Convert sampling vector to integer to extract correct elements from simulation
        samps1 = convert.(Int, t2cor);

        try
            
            mm11,mm22 = 0, 0
            try
                # Simulate
                SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllCp_Scall(ts1, pD2, ivss1, samps1);
                
                mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat["Erros"][:,i,1].^2) .+ ((SimOnTime1[:,14] .- dat["Means"][:,i,1]).^2)./(dat["Erros"][:,i,1].^2)));
                mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat["Erros"][:,i,2].^2) .+ (((SimOnTime1[:,13]) .- dat["Means"][:,i,2]).^2)./(dat["Erros"][:,i,2].^2)));

                
            catch

                # Simulate
                t2cor = dat["sts"][1:end-1,i]
                ts1 = collect(0:t2cor[end]);
                samps1 = convert.(Int, t2cor);
                SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllCp_Scall(ts1, pD2, ivss1, samps1);

                mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat["Erros"][1:end-1,i,1].^2) .+ ((SimOnTime1[:,14] .- dat["Means"][1:end-1,i,1]).^2)./(dat["Erros"][1:end-1,i,1].^2)));
                mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat["Erros"][1:end-1,i,2].^2) .+ (((SimOnTime1[:,13]) .- dat["Means"][1:end-1,i,2]).^2)./(dat["Erros"][1:end-1,i,2].^2)));
            end
                
            obj = -(mm22)+sqrt(abs(mm11))+obj;
            
        catch
            obj=1e30
        end

        
    end
    
    return(obj)

end



################################################################################################

function genStanInitDict2(samps, names)
    chains = size(repeat(preTrans, 1, 10))[2]
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


###############################################################################################



function CompRepUtil_Entro_OEDmc2(ins)

    # Input converted in the write format
    Pyr = round(ins[1], digits = 1)*1000;

    # Simulate the model for each theta and Y0 sampled
    SimulsAll = Array{Any}(undef, 1000, 3);
    for i in 1:1000
        ivss = [Pyr, 0, 0, sampsY0[i,1], sampsY0[i,2], sampsY0[i,3], 0, 0, 0, 0, 0, Pyr, 0, Pyr*(sampsTh[1,end]*0.12)];
        SimulsAll[i,1], SimulsAll[i,2], SimulsAll[i,3]  = PyruvateHP_NMR_SolveAllCp_Scall(ts, vcat(sampsTh[i,:], 0), ivss, samps);
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
            fitts1[names[1]] = Distributions.fit.(dists, Ref(Obs[j,:]));
            bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2];
            bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]];
        catch
            fitts1[names[1]] = Distributions.fit.(dists2, Ref(Obs[j,:]));
            bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]+1;
            bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs[j,:])))[2]];
        end

        EntObs[j] = entropy(bestfit1[names[1]]); 
    end

    HES = zeros(1,1);
    HES[1,1] = sum(EntObs)

    return(HES[1,1])
end


######################################################################################


function simDwnExp2(expS, datsD, PyrCons, nSamps, poster, allD, plt = false)
    tsC2 = vcat(collect(0:5:395));
    ts = 0:tsC2[end];
    samps = convert.(Int, tsC2);


    exin = findall(expS .== datsD)
    postDwn = poster[StatsBase.sample(1:size(poster)[1], nSamps, replace = false),:]

    scF = allD[exin][1][1,2]/(PyrCons[1]*1000*0.12)
    simsPrior = Array{Any}(undef, nSamps);
    for i in 1:nSamps
        ivss = [round(PyrCons[1], digits = 1)*1000, 0, 0, postDwn[i,15], postDwn[i,16], postDwn[i,17], 0, 0, 0, 0, 0, round(PyrCons[1], digits = 1)*1000, 0, round(PyrCons[1], digits = 1)*1000*(scF*0.12)];;
        samps = convert.(Int, tsC2);

        SimulsAll = []
        
        SimulsAll, SimulsAll2, SimulsAll3 = PyruvateHP_NMR_SolveAllCp_Scall(ts, vcat(postDwn[i,1:13], scF, 0), ivss, samps);

        simsPrior[i] = SimulsAll;
    end

    # Dims = Up/Down, Observable, Experiment
    PriorQuant = Array{Any}(undef,2,1,1);
    lacPirorSims = zeros(size(simsPrior[1])[1], nSamps)
    [lacPirorSims[:,i] = simsPrior[i][:,13,1] for i in 1:nSamps ];
    for k in 1:1
        for m in 1:1
            PriorQuant[1,k,m] = [percentile(lacPirorSims[j,:], 99.5) for j in 1:size(lacPirorSims)[1]]; # Up
            PriorQuant[2,k,m] = [percentile(lacPirorSims[j,:], 0.5 ) for j in 1:size(lacPirorSims)[1]]; # Down
        end
    end

    if plt == true
        pr = scatter(allD[exin][1][:,1]./60, allD[exin][1][:,4], yerror = allD[exin][1][:,5], label = "")
        plot!(tsC2./60, PriorQuant[1], label = "", color = "red", linewidth = 3,
            margin=10Plots.mm,xtickfont=font(16), ytickfont=font(16), guidefont=font(18), titlefont=font(18),
            xlabel = "time (min)", ylabel = "Lactate (AU)")
        plot!(tsC2./60, PriorQuant[2], label = "", color = "red", linewidth = 3)
        plot!(tsC2./60, PriorQuant[1], fillrange=PriorQuant[2], label="", color="red", fillalpha=0.2)

        display(pr)
    end

    return(simsPrior, PriorQuant)

end





######################################################################################







######################################################################################







######################################################################################







######################################################################################







######################################################################################







######################################################################################







######################################################################################







######################################################################################