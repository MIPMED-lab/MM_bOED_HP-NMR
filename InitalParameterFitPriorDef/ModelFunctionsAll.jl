

# -------------------------------------------------------------------- NO TRANSPORTER --------------------------------------------------------------------

function PyruvateHP_Cells!(du, u, p, t)
    # Php, Xhp, P, X = u;
    # k_f1, k_r1, T1_P, T1_X = p;

    Php, Xhp = u
    k_f1, k_r1, T1_X, T1_P = p

    du[1] = dPhp = k_r1 * Xhp - ((1 / T1_P) + k_f1) * Php
    du[2] = dXhp = k_f1 * Php - ((1 / T1_X) + k_r1) * Xhp

end



function PyruvateHP_NMR_SolveAll(ts, pD, ivss, samps)

    if length(size(pD)) == 1
        pD = reshape(pD, size(pD)[1], 1)
    end

    if size(pD)[2] != 5
        pD = pD'
    end

    if length(ivss) / 2 > 1
        if size(ivss)[2] != 2
            ivss = ivss'
        end
    end

    AllSolTest = zeros(length(samps), 2, length(pD[:, 1])) # Simulation of the system observed
    AllSolTest_Off = Array{Any,1}(undef, length(pD[:, 1])) # Simulation of the system before we obvserve it (considering time offset). First column is the time vector
    AllSolTest_Tog = Array{Any,1}(undef, length(pD[:, 1])) # Previous two together. First column is the time vector. 

    for drawInd in collect(1:length(pD[:, 1]))

        p = pD[drawInd, 1:end-1]
        tau = pD[drawInd, end]

        if length(ivss) / 2 > 1
            ivss2 = ivss[drawInd, :]
        else
            ivss2 = ivss
        end

        # Account for the time offset

        prob_off = ODEProblem(PyruvateHP_Cells!, ivss2, (-tau, 0), p)
        part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(), reltol=1.0e-9, abstol=1.0e-9)

        ivss2 = part1_off.u[end]

        prob = ODEProblem(PyruvateHP_Cells!, ivss2, (ts[1], ts[end]), p)
        part1 = DifferentialEquations.solve(prob, CVODE_BDF(), reltol=1.0e-9, abstol=1.0e-9, saveat=1)
        

        AllSolTest[:, :, drawInd]

        tmp = zeros(length(part1_off.u), 2)



        for j in 1:2
            AllSolTest[:, j, drawInd] = [part1.u[i][j] for i in 1:length(part1.u)][samps.+1]
            tmp[:, j] = [part1_off.u[i][j] for i in 1:length(part1_off.u)]
        end

        AllSolTest_Off[drawInd] = hcat(part1_off.t, tmp)
        AllSolTest_Tog[drawInd] = vcat(AllSolTest_Off[drawInd][1:end-1, :], hcat(samps, AllSolTest[:, :, drawInd]))
        




    end

    return AllSolTest, AllSolTest_Off, AllSolTest_Tog

end

# Objective Function No Transporter

function ObjectFunctME(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    pD1 = vcat(vcat(p[1:2].*2, p[3:end]), 0);
    pD2 = vcat(vcat(p[1:2].*4, p[3:end]), 0);
    pD3 = vcat(vcat(p[1:2].*8, p[3:end]), 0);

    # Define time vector
    t2cor = dat2[:,1];
    t4cor = dat4[:,1];
    t8cor = dat8[:,1];

    # Define equaly-spaced time vector
    ts1 = collect(0:t2cor[end]);
    ts2 = collect(0:t4cor[end]);
    ts3 = collect(0:t8cor[end]);

    # Define initial value for simulation (use of experimental mean)
    ivss1 = [dat2[1,2], dat2[1,4]];
    ivss2 = [dat4[1,2], dat4[1,4]];
    ivss3 = [dat8[1,2], dat8[1,4]];

    # Convert sampling vector to integer to extract correct elements from simulation
    samps1 = convert.(Int, t2cor);
    samps2 = convert.(Int, t4cor);
    samps3 = convert.(Int, t8cor);


    # Simulate
    SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAll(ts1, pD1, ivss1, samps1);
    SimOnTime2, SimOffTime2, SimAll2  = PyruvateHP_NMR_SolveAll(ts2, pD2, ivss2, samps2);
    SimOnTime3, SimOffTime3, SimAll3  = PyruvateHP_NMR_SolveAll(ts3, pD3, ivss3, samps3);


    # Use of log-likelihood as cost funtion: J_llk = sum(-1/2 * (log(2*pi) + log(std^2) + (sim-dat)^2/std^2))
    mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat2[:,3].^2) .+ ((SimOnTime1[:,1] .- dat2[:,2]).^2)./(dat2[:,3].^2)));
    mm12 = sum((-1/2) .* (log(2*pi) .+ log.(dat4[:,3].^2) .+ ((SimOnTime2[:,1] .- dat4[:,2]).^2)./(dat4[:,3].^2)));
    mm13 = sum((-1/2) .* (log(2*pi) .+ log.(dat8[:,3].^2) .+ ((SimOnTime3[:,1] .- dat8[:,2]).^2)./(dat8[:,3].^2)));

    mm21 = sum((-1/2) .* (log(2*pi) .+ log.(dat2[:,5].^2) .+ ((SimOnTime1[:,2] .- dat2[:,4]).^2)./(dat2[:,5].^2)));
    mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat4[:,5].^2) .+ ((SimOnTime2[:,2] .- dat4[:,4]).^2)./(dat4[:,5].^2)));
    mm23 = sum((-1/2) .* (log(2*pi) .+ log.(dat8[:,5].^2) .+ ((SimOnTime3[:,2] .- dat8[:,4]).^2)./(dat8[:,5].^2)));


    obj = (mm21+mm22+mm23)*(mm11+mm12+mm13);
        
    return(obj)

end

# -------------------------------------------------------------------- TRANSPORTER --------------------------------------------------------------------

function PyruvateHP_CellsTb!(du, u, p, t)

    Pout, Php, Xhp = u;
    k_f1, k_r1, T1_X, T1_P, kin = p;


    du[1] = dPout = -Pout*kin - (Pout/T1_P)
    du[2] = dPhp = k_r1*Xhp - ((1/T1_P)+k_f1)*Php + Pout*kin;
    du[3] = dXhp = k_f1*Php - ((1/T1_X)+k_r1)*Xhp;

end


function PyruvateHP_NMR_SolveAllTb(ts, pD, ivss, samps)

    if length(size(pD)) == 1
        pD = reshape(pD,size(pD)[1],1);
    end

    if size(pD)[2] != 6     
        pD = pD';
    end

    if length(ivss)/3 > 1
        if size(ivss)[2] != 3
            ivss = ivss';
        end
    end

    AllSolTest = zeros(length(samps), 3, length(pD[:,1])); # Simulation of the system observed
    AllSolTest_Off = Array{Any,1}(undef,length(pD[:,1])); # Simulation of the system before we obvserve it (considering time offset). First column is the time vector
    AllSolTest_Tog = Array{Any,1}(undef,length(pD[:,1])); # Previous two together. First column is the time vector. 


    
    for drawInd in collect(1:length(pD[:,1]))
        
        p = pD[drawInd,1:end-1];
        tau = pD[drawInd,end];
        
        if length(ivss)/3 > 1
            ivss2 = ivss[drawInd,:];
        else
            ivss2 = ivss;
        end
    
        # Account for the time offset
    
        
            prob_off = ODEProblem(PyruvateHP_CellsTb!,ivss2,(-tau, 0),p);
            part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9);
        
            ivss2 = part1_off.u[end];
        
            prob = ODEProblem(PyruvateHP_CellsTb!,ivss2,(ts[1], ts[end]),p);
            part1 = DifferentialEquations.solve(prob, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
        
        
        AllSolTest[:,:,drawInd]
    
        tmp = zeros(length(part1_off.u), 3);
    
        
        

        for j in 1:3
            AllSolTest[:,j,drawInd] = [part1.u[i][j] for i in 1:length(part1.u)][samps.+1];
            tmp[:,j] = [part1_off.u[i][j] for i in 1:length(part1_off.u)];
        end

        AllSolTest_Off[drawInd] = hcat(part1_off.t, tmp);
        AllSolTest_Tog[drawInd] = vcat(AllSolTest_Off[drawInd][1:end-1, :], hcat(samps, AllSolTest[:,:,drawInd]));
        
        
        
    
    end

    return AllSolTest, AllSolTest_Off, AllSolTest_Tog

end



function ObjectFunctMETR(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    pD1 = vcat(vcat(p[1:2].*2, p[3:end-1], p[end]*2), 0);
    pD2 = vcat(vcat(p[1:2].*4, p[3:end-1], p[end]*4), 0);
    pD3 = vcat(vcat(p[1:2].*8, p[3:end-1], p[end]*8), 0);

    # Define time vector
    t2cor = dat2[:,1];
    t4cor = dat4[:,1];
    t8cor = dat8[:,1];

    # Define equaly-spaced time vector
    ts1 = collect(0:t2cor[end]);
    ts2 = collect(0:t4cor[end]);
    ts3 = collect(0:t8cor[end]);

    # Define initial value for simulation (use of experimental mean)
    ivss1 = [dat2[1,2], 0, dat2[1,4]];
    ivss2 = [dat4[1,2], 0, dat4[1,4]];
    ivss3 = [dat8[1,2], 0, dat8[1,4]];

    # Convert sampling vector to integer to extract correct elements from simulation
    samps1 = convert.(Int, t2cor);
    samps2 = convert.(Int, t4cor);
    samps3 = convert.(Int, t8cor);


    # Simulate
    SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllTb(ts1, pD1, ivss1, samps1);
    SimOnTime2, SimOffTime2, SimAll2  = PyruvateHP_NMR_SolveAllTb(ts2, pD2, ivss2, samps2);
    SimOnTime3, SimOffTime3, SimAll3  = PyruvateHP_NMR_SolveAllTb(ts3, pD3, ivss3, samps3);


    # Use of log-likelihood as cost funtion: J_llk = sum(-1/2 * (log(2*pi) + log(std^2) + (sim-dat)^2/std^2))
    mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat2[:,3].^2) .+ (((SimOnTime1[:,1]+SimOnTime1[:,2]) .- dat2[:,2]).^2)./(dat2[:,3].^2)));
    mm12 = sum((-1/2) .* (log(2*pi) .+ log.(dat4[:,3].^2) .+ (((SimOnTime2[:,1]+SimOnTime2[:,2]) .- dat4[:,2]).^2)./(dat4[:,3].^2)));
    mm13 = sum((-1/2) .* (log(2*pi) .+ log.(dat8[:,3].^2) .+ (((SimOnTime3[:,1]+SimOnTime3[:,2]) .- dat8[:,2]).^2)./(dat8[:,3].^2)));

    mm21 = sum((-1/2) .* (log(2*pi) .+ log.(dat2[:,5].^2) .+ ((SimOnTime1[:,3] .- dat2[:,4]).^2)./(dat2[:,5].^2)));
    mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat4[:,5].^2) .+ ((SimOnTime2[:,3] .- dat4[:,4]).^2)./(dat4[:,5].^2)));
    mm23 = sum((-1/2) .* (log(2*pi) .+ log.(dat8[:,5].^2) .+ ((SimOnTime3[:,3] .- dat8[:,4]).^2)./(dat8[:,5].^2)));


    obj = (mm21+mm22+mm23)*(mm11+mm12+mm13);
        
    return(obj)

end






# -------------------------------------------------------------------- REPRESSION --------------------------------------------------------------------

function PyruvateHP_CellsRp!(du, u, p, t)

    Pout, Php, Xhp, LDH, LDHna = u;
    k_f1, k_r1, T1_X, T1_P, kin, kf, kr = p;


    du[1] = dPout = - (Pout*kin) - (Pout/T1_P)
    du[2] = dPhp = (Pout*kin) + (k_r1*Xhp*LDH) + (kr*LDHna) - (k_f1*Php*LDH) - (kf*LDH*Php) - (Php/T1_P)
    du[3] = dXhp = (k_f1*Php*LDH) - (k_r1*Xhp*LDH) - (Xhp/T1_X)
    # du[2] = dPhp = k_r1*Xhp*LDH + kr * LDHna - ((1/T1_P)+k_f1)*Php*LDH - (kf * LDH * Php) + Pout*kin;
    # du[3] = dXhp = k_f1*Php*LDH - ((1/T1_X)+k_r1)*Xhp*LDH;

    du[4] = dLDH   = (kr * LDHna) - (kf * LDH * Php) 
    du[5] = dLDHna = (kf * LDH * Php) - (kr * LDHna)

end


function PyruvateHP_NMR_SolveAllRp(ts, pD, ivss, samps)

    if length(size(pD)) == 1
        pD = reshape(pD,size(pD)[1],1);
    end

    if size(pD)[2] != 8     
        pD = pD';
    end

    if length(ivss)/5 > 1
        if size(ivss)[2] != 5
            ivss = ivss';
        end
    end

    AllSolTest = zeros(length(samps), 5, length(pD[:,1])); # Simulation of the system observed
    AllSolTest_Off = Array{Any,1}(undef,length(pD[:,1])); # Simulation of the system before we obvserve it (considering time offset). First column is the time vector
    AllSolTest_Tog = Array{Any,1}(undef,length(pD[:,1])); # Previous two together. First column is the time vector. 


    
    for drawInd in collect(1:length(pD[:,1]))
        
        p = pD[drawInd,1:end-1];
        tau = pD[drawInd,end];
        
        if length(ivss)/5 > 1
            ivss2 = ivss[drawInd,:];
        else
            ivss2 = ivss;
        end
    
        # Account for the time offset
    
        
            prob_off = ODEProblem(PyruvateHP_CellsRp!,ivss2,(-tau, 0),p);
            part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9);
        
            ivss2 = part1_off.u[end];
        
            prob = ODEProblem(PyruvateHP_CellsRp!,ivss2,(ts[1], ts[end]),p);
            part1 = DifferentialEquations.solve(prob, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
        
        
        AllSolTest[:,:,drawInd]
    
        tmp = zeros(length(part1_off.u), 5);
    
        
        

        for j in 1:5
            AllSolTest[:,j,drawInd] = [part1.u[i][j] for i in 1:length(part1.u)][samps.+1];
            tmp[:,j] = [part1_off.u[i][j] for i in 1:length(part1_off.u)];
        end

        AllSolTest_Off[drawInd] = hcat(part1_off.t, tmp);
        AllSolTest_Tog[drawInd] = vcat(AllSolTest_Off[drawInd][1:end-1, :], hcat(samps, AllSolTest[:,:,drawInd]));
        
        
        
    
    end

    return AllSolTest, AllSolTest_Off, AllSolTest_Tog

end


function ObjectFunctMERP(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    pD2 = vcat(vcat(p[1:2].*4, p[3:4], p[5:end-1]*4), 0);

    # Define time vector
    t2cor = dat32mM[:,1];
    t4cor = dat128mM[:,1];

    # Define equaly-spaced time vector
    ts1 = collect(0:t2cor[end]);
    ts2 = collect(0:t4cor[end]);

    # Define initial value for simulation (use of experimental mean)
    ivss1 = [dat32mM[1,2], 0, dat32mM[1,4], p[end]*4, 0];
    ivss2 = [dat128mM[1,2], 0, dat128mM[1,4], p[end]*4, 0];

    # Convert sampling vector to integer to extract correct elements from simulation
    samps1 = convert.(Int, t2cor);
    samps2 = convert.(Int, t4cor);

    # Simulate
    SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllRp(ts1, pD2, ivss1, samps1);
    SimOnTime2, SimOffTime2, SimAll2  = PyruvateHP_NMR_SolveAllRp(ts2, pD2, ivss2, samps2);

    # Use of log-likelihood as cost funtion: J_llk = sum(-1/2 * (log(2*pi) + log(std^2) + (sim-dat)^2/std^2))
    mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,3].^2) .+ (((SimOnTime1[:,1]+SimOnTime1[:,2]) .- dat32mM[:,2]).^2)./(dat32mM[:,3].^2)));
    mm12 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,3].^2) .+ (((SimOnTime2[:,1]+SimOnTime2[:,2]) .- dat128mM[:,2]).^2)./(dat128mM[:,3].^2)));

    mm21 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,5].^2) .+ ((SimOnTime1[:,3] .- dat32mM[:,4]).^2)./(dat32mM[:,5].^2)));
    mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,5].^2) .+ ((SimOnTime2[:,3] .- dat128mM[:,4]).^2)./(dat128mM[:,5].^2)));


    obj = (mm21+mm22)*(mm11+mm12);
        
    return(obj)

end





# -------------------------------------------------------------------- ALLOSTERIC --------------------------------------------------------------------

function PyruvateHP_CellsAl!(du, u, p, t)

    Pout, Php, Xhp, LDH, LDHP1, LDHP2, LDHP3, LDHP4, XhpHyper, PhpHyper, Obs_XhpHyper, Obs_PhpHyper = u;
    T1_X, T1_P, kin, kpl, kf, kr, a1, a2, a3, b1, b2, b3, ScFm = p;


    du[1] = dPout = - (Pout*kin)
    du[2] = dPhp = (Pout*kin) + (kr*LDHP1) + (a1*kr*LDHP2) + (a2*kr*LDHP3) + (a3*kr*LDHP4) - (kf*LDH*Php) - (a1*kf*LDHP1*Php) - (a2*kf*LDHP2*Php) - (a3*kf*LDHP3*Php) 
    du[3] = dXhp = (kpl*LDHP1)*3 + (b1*kpl*LDHP2)*2 + (b2*kpl*LDHP3)*2 + (b3*kpl*LDHP4)
    du[4] = dLDH = (kr*LDHP1) + (kpl*LDHP1)*3 - (kf*LDH*Php)
    du[5] = dLDHP1 = (kf*LDH*Php) + (a1*kr*LDHP2) + (b1*kpl*LDHP2)*2 - (kr*LDHP1) - (a1*kf*LDHP1*Php) - (kpl*LDHP1)*3
    du[6] = dLDHP2 = (a1*kf*LDHP1*Php) + (a2*kr*LDHP3) + (b2*kpl*LDHP3)*2 - (a1*kr*LDHP2) - (a2*kf*LDHP2*Php) - (b1*kpl*LDHP2)*2
    du[7] = dLDHP3 = (a2*kf*LDHP2*Php) + (a3*kr*LDHP4) + (b3*kpl*LDHP4) - (a2*kr*LDHP3) - (a3*kf*LDHP3*Php) - (b2*kpl*LDHP3)*2
    du[8] = dLDHP4 = (a3*kf*LDHP3*Php) - (a3*kr*LDHP4) - (b3*kpl*LDHP4)

    du[9] = dXhpHyper = dXhp - (XhpHyper/T1_X)
    du[10] = dPhpHyper = dPout + dPhp  + dLDHP1 + dLDHP2*2 + dLDHP3*3 + dLDHP4*4  - (PhpHyper/T1_P)

    du[11] = dObs_XhpHyper =  dXhpHyper*(ScFm*0.12)
    du[12] = dObs_PhpHyper =  dPhpHyper*(ScFm*0.12)


end



function PyruvateHP_NMR_SolveAllAl(ts, pD, ivss, samps)


    nStat = 12;
    nPar = 14;

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
    
        
            prob_off = ODEProblem(PyruvateHP_CellsAl!,ivss2,(-tau, 0),p);
            part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(),reltol=1.0e-14,abstol=1.0e-14);
        
            ivss2 = part1_off.u[end];
        
            prob = ODEProblem(PyruvateHP_CellsAl!,ivss2,(ts[1], ts[end]),p);
            part1 = DifferentialEquations.solve(prob, CVODE_BDF(),reltol=1.0e-14,abstol=1.0e-14,saveat=1);
        
        
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



function ObjectFunctMEAl(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    pD2 = vcat(p[1:end-1], ScFm, 0);

    # Define time vector
    t2cor = dat32mM[:,1];
    t4cor = dat128mM[:,1];

    # Define equaly-spaced time vector
    ts1 = collect(0:t2cor[end]);
    ts2 = collect(0:t4cor[end]);

    # Define initial value for simulation (use of experimental mean)
    ivss1 = [dat32mM[1,2]/(ScFm*0.12), 0, dat32mM[1,4]/(ScFm*0.12), p[end], 0, 0, 0, 0, 0, dat32mM[1,2]/(ScFm*0.12), dat32mM[1,4], dat32mM[1,2]];
    ivss2 = [dat128mM[1,2]/(ScFm*0.12), 0, dat128mM[1,4]/(ScFm*0.12), p[end], 0, 0, 0, 0, 0, dat128mM[1,2]/(ScFm*0.12), dat128mM[1,4], dat128mM[1,2]];

    # Convert sampling vector to integer to extract correct elements from simulation
    samps1 = convert.(Int, t2cor);
    samps2 = convert.(Int, t4cor);

    # Simulate
    SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllAl(ts1, pD2, ivss1, samps1);
    SimOnTime2, SimOffTime2, SimAll2  = PyruvateHP_NMR_SolveAllAl(ts2, pD2, ivss2, samps2);

    # Use of log-likelihood as cost funtion: J_llk = sum(-1/2 * (log(2*pi) + log(std^2) + (sim-dat)^2/std^2))
    # mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,3].^2) .+ (((SimOnTime1[:,1]+SimOnTime1[:,2]) .- dat32mM[:,2]).^2)./(dat32mM[:,3].^2)));
    # mm12 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,3].^2) .+ (((SimOnTime2[:,1]+SimOnTime2[:,2]) .- dat128mM[:,2]).^2)./(dat128mM[:,3].^2)));

    mm21 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,5].^2) .+ ((SimOnTime1[:,11] .- dat32mM[:,4]).^2)./(dat32mM[:,5].^2)));
    mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,5].^2) .+ ((SimOnTime2[:,11] .- dat128mM[:,4]).^2)./(dat128mM[:,5].^2)));


    # obj = (mm21+mm22)*(mm11+mm12);
    obj = -(mm21+mm22);

    
    return(obj)

end




# -------------------------------------------------------------------- COMPETITIVE REPRESSION --------------------------------------------------------------------

function PyruvateHP_CellsCp!(du, u, p, t)

    Pout, Php, Xhp, NADH, NAD, LDH, LDHac, LDHacE, LDHna, LDHP, XhpHyper, PhpHyper, Obs_XhpHyper, Obs_PhpHyper = u;
    T1_X, T1_P, kin, kpl, kbn, kun, kbp, kup, kunE, kbnE, ki, kr, ScFm = p;


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
    du[11] = dXhpHyper =  dXhp - (XhpHyper/T1_X)  #  These are transition states to compute the next two ones (A.U. too) assuming that all is hyperpolarised. Tecchnically it could be coupled with the next two ones for a final result too. 
    du[12] = dPhpHyper =  dPout + dPhp + dLDHna + dLDHP - (PhpHyper/T1_P)

    du[13] = dObs_XhpHyper =  dXhpHyper * (ScFm*0.12)
    du[14] = dObs_PhpHyper =  dPhpHyper * (ScFm*0.12)

end


function PyruvateHP_NMR_SolveAllCp(ts, pD, ivss, samps)


    nStat = 14;
    nPar = 14;

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
    
        
            prob_off = ODEProblem(PyruvateHP_CellsCp!,ivss2,(-tau, 0),p);
            part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9);
        
            ivss2 = part1_off.u[end];
        
            prob = ODEProblem(PyruvateHP_CellsCp!,ivss2,(ts[1], ts[end]),p);
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



function ObjectFunctMECp(p)

    # Define parameter vectors for each amount of cells (last parameter is a time delay, not used in here)
    pD2 = vcat(p[1:end-3], ScFm, 0);

    # Define time vector
    t2cor = dat32mM[:,1];
    t4cor = dat128mM[:,1];

    # Define equaly-spaced time vector
    ts1 = collect(0:t2cor[end]);
    ts2 = collect(0:t4cor[end]);

    # Define initial value for simulation (use of experimental mean)
    ivss1 = [dat32mM[1,2]/(ScFm*0.12), 0, dat32mM[1,4]/(ScFm*0.12), p[end-2], p[end-1], p[end], 0, 0, 0, 0, 0, dat32mM[1,2]/(ScFm*0.12), dat32mM[1,4], dat32mM[1,2]];
    ivss2 = [dat128mM[1,2]/(ScFm*0.12), 0, dat128mM[1,4]/(ScFm*0.12), p[end-2], p[end-1], p[end], 0, 0, 0, 0, 0, dat128mM[1,2]/(ScFm*0.12), dat128mM[1,4], dat128mM[1,2]];

    # Convert sampling vector to integer to extract correct elements from simulation
    samps1 = convert.(Int, t2cor);
    samps2 = convert.(Int, t4cor);

    # Simulate
    SimOnTime1, SimOffTime1, SimAll1  = PyruvateHP_NMR_SolveAllCp(ts1, pD2, ivss1, samps1);
    SimOnTime2, SimOffTime2, SimAll2  = PyruvateHP_NMR_SolveAllCp(ts2, pD2, ivss2, samps2);

    # Use of log-likelihood as cost funtion: J_llk = sum(-1/2 * (log(2*pi) + log(std^2) + (sim-dat)^2/std^2))
    # mm11 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,3].^2) .+ (((SimOnTime1[:,1]+SimOnTime1[:,2]) .- dat32mM[:,2]).^2)./(dat32mM[:,3].^2)));
    # mm12 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,3].^2) .+ (((SimOnTime2[:,1]+SimOnTime2[:,2]) .- dat128mM[:,2]).^2)./(dat128mM[:,3].^2)));

    # mm21 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,5].^2) .+ ((SimOnTime1[:,11] .- dat32mM[:,4]).^2)./(dat32mM[:,5].^2)));
    # mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,5].^2) .+ ((SimOnTime2[:,11] .- dat128mM[:,4]).^2)./(dat128mM[:,5].^2)));

    # mm21 = sum((-1/2) .* (log(2*pi) .+ log.(1^2) .+ ((SimOnTime1[:,13] .- dat32mM[:,4]).^2)./(1^2)));
    # mm22 = sum((-1/2) .* (log(2*pi) .+ log.(1^2) .+ ((SimOnTime2[:,13] .- dat128mM[:,4]).^2)./(1^2)));
    
    mm21 = sum((-1/2) .* (log(2*pi) .+ log.(dat32mM[:,5].^2) .+ ((SimOnTime1[:,13] .- dat32mM[:,4]).^2)./(dat32mM[:,5].^2)));
    mm22 = sum((-1/2) .* (log(2*pi) .+ log.(dat128mM[:,5].^2) .+ ((SimOnTime2[:,13] .- dat128mM[:,4]).^2)./(dat128mM[:,5].^2)));

    # obj = (mm21+mm22)*(mm11+mm12);
    obj = -(mm21+mm22);

    
    return(obj)

end




# -------------------------------------------------------------------- Michaelis Menten Kinetics --------------------------------------------------------------------
function PyruvateHP_CellsMM!(du, u, p, t)

    Pout, Pin, LDH, LDHPyr, Lac = u;
    kin, kf, kr, kpl, T1_P, T1_L = p;

    du[1] = dPout = -Pout*kin - Pout/T1_P
    du[2] = dPin = Pout*kin + kr*LDHPyr - kf*Pin*LDH - Pin/T1_P
    du[3] = dLDH = kr*LDHPyr + kpl*LDHPyr - kf*Pin*LDH
    du[4] = dLDHPyr = kf*Pin*LDH - kr*LDHPyr - kpl*LDHPyr - LDHPyr/T1_P
    du[5] = dLac = kpl*LDHPyr - (Lac/T1_L)

end


function PyruvateHP_NMR_SolveAllMM(ts, pD, ivss, samps)

    if length(size(pD)) == 1
        pD = reshape(pD,size(pD)[1],1);
    end

    if size(pD)[2] != 7     
        pD = pD';
    end

    if length(ivss)/5 > 1
        if size(ivss)[2] != 5
            ivss = ivss';
        end
    end

    AllSolTest = zeros(length(samps), 5, length(pD[:,1])); # Simulation of the system observed
    AllSolTest_Off = Array{Any,1}(undef,length(pD[:,1])); # Simulation of the system before we obvserve it (considering time offset). First column is the time vector
    AllSolTest_Tog = Array{Any,1}(undef,length(pD[:,1])); # Previous two together. First column is the time vector. 


    
    for drawInd in collect(1:length(pD[:,1]))
        
        p = pD[drawInd,1:end-1];
        tau = pD[drawInd,end];
        
        if length(ivss)/5 > 1
            ivss2 = ivss[drawInd,:];
        else
            ivss2 = ivss;
        end
    
        # Account for the time offset
    
            prob_off = ODEProblem(PyruvateHP_CellsMM!,ivss2,(-tau, 0),p);
            part1_off = DifferentialEquations.solve(prob_off, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9);
        
            ivss2 = part1_off.u[end];
        
            prob = ODEProblem(PyruvateHP_CellsMM!,ivss2,(ts[1], ts[end]),p);
            part1 = DifferentialEquations.solve(prob, CVODE_BDF(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
        
        
        AllSolTest[:,:,drawInd]
    
        tmp = zeros(length(part1_off.u), 5);

        for j in 1:5
            AllSolTest[:,j,drawInd] = [part1.u[i][j] for i in 1:length(part1.u)][samps.+1];
            tmp[:,j] = [part1_off.u[i][j] for i in 1:length(part1_off.u)];
        end

        AllSolTest_Off[drawInd] = hcat(part1_off.t, tmp);
        AllSolTest_Tog[drawInd] = vcat(AllSolTest_Off[drawInd][1:end-1, :], hcat(samps, AllSolTest[:,:,drawInd]));
    
    end

    return AllSolTest, AllSolTest_Off, AllSolTest_Tog

end