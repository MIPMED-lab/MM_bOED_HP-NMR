function restructDatInf(paths)

    ######################## OBSERVABLES ########################

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

    Y0us = Array{Float64,2}(undef, 3, m);
    tml = maximum([length(DataAll[i][1,1]:1:DataAll[i][end,1]) for i in 1:m]);
    ts = Array{Float64,2}(undef, tml, m);
    itp = Array{Float64,1}(undef, m);

    ms = collect(2:2:obser*2);
    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        Y0us[:,i] = [DataAll[i][1,2], 0, 0]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function
        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
        if occursin("2Mil", paths[i])
            ncells[1,i] = 2;
        elseif occursin("4Mil", paths[i])
            ncells[1,i] = 4;
        elseif occursin("8Mil", paths[i])
            ncells[1,i] = 8;
        end
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
        "nts" => nts
    
    );

    return data_multi

end






function restructDatInfNT(paths)

    ######################## OBSERVABLES ########################

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

    Y0us = Array{Float64,2}(undef, 2, m);
    tml = maximum([length(DataAll[i][1,1]:1:DataAll[i][end,1]) for i in 1:m]);
    ts = Array{Float64,2}(undef, tml, m);
    itp = Array{Float64,1}(undef, m);

    ms = collect(2:2:obser*2);
    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        Y0us[:,i] = [DataAll[i][1,2], 0]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function
        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
        if occursin("2Mil", paths[i])
            ncells[1,i] = 2;
        elseif occursin("4Mil", paths[i])
            ncells[1,i] = 4;
        elseif occursin("8Mil", paths[i])
            ncells[1,i] = 8;
        end
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
        "nts" => nts
    
    );

    return data_multi

end



function restructDatInfAllo(paths, LDHmean)

    ######################## OBSERVABLES ########################

    P0s = JLD2.load("C:\\IBECPostDocDrive\\2024_01_16_NCvsKR\\DataProcessingInference\\ProcessedData\\P0s.jld")["P0s"][1:11];
    ScFm = mean(P0s./(3200*0.12));

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

    Y0us = Array{Float64,2}(undef, 12, m);
    tml = maximum([length(DataAll[i][1,1]:1:DataAll[i][end,1]) for i in 1:m]);
    ts = Array{Float64,2}(undef, tml, m);
    itp = Array{Float64,1}(undef, m);

    ms = collect(2:2:obser*2);
    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        Y0us[:,i] = [DataAll[i][1,2]/(ScFm*0.12), 0, 0, LDHmean, 0, 0, 0, 0, 0, DataAll[i][1,2]/(ScFm*0.12), 0, DataAll[i][1,2]]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
        if occursin("2Mil", paths[i])
            ncells[1,i] = 2;
        elseif occursin("4Mil", paths[i])
            ncells[1,i] = 4;
        elseif occursin("8Mil", paths[i])
            ncells[1,i] = 8;
        end
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
        "nts" => nts
    
    );

    return data_multi

end




function restructDatInfCompRep(paths, NADHmean, NADmean, LDHmean)

    ######################## OBSERVABLES ########################

    P0s = JLD2.load("C:\\IBECPostDocDrive\\2024_01_16_NCvsKR\\DataProcessingInference\\ProcessedData\\P0s.jld")["P0s"][1:11];
    ScFm = mean(P0s./(3200*0.12));

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
    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        Y0us[:,i] = [DataAll[i][1,2]/(ScFm*0.12), 0,0, NADHmean, NADmean, LDHmean, 0, 0, 0, 0, 0, DataAll[i][1,2]/(ScFm*0.12), 0, DataAll[i][1,2]]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
        if occursin("2Mil", paths[i])
            ncells[1,i] = 2;
        elseif occursin("4Mil", paths[i])
            ncells[1,i] = 4;
        elseif occursin("8Mil", paths[i])
            ncells[1,i] = 8;
        end
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
        "nts" => nts
    
    );

    return data_multi

end




function restructDatInfCompRep_FixedSF(paths, NADHmean, NADmean, LDHmean, ScFm)

    ######################## OBSERVABLES ########################

    # P0s = JLD2.load("C:\\IBECPostDocDrive\\2024_01_16_NCvsKR\\DataProcessingInference\\ProcessedData\\P0s.jld")["P0s"][1:11];
    # ScFm = mean(P0s./(3200*0.12));

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
    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        Y0us[:,i] = [DataAll[i][1,2]/(ScFm*0.12), 0,0, NADHmean, NADmean, LDHmean, 0, 0, 0, 0, 0, DataAll[i][1,2]/(ScFm*0.12), 0, DataAll[i][1,2]]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
        if occursin("2Mil", paths[i])
            ncells[1,i] = 2;
        elseif occursin("4Mil", paths[i])
            ncells[1,i] = 4;
        elseif occursin("8Mil", paths[i])
            ncells[1,i] = 8;
        end
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
        "nts" => nts
    
    );

    return data_multi

end




function restructDatInfCompRep_RealExpDat(paths, conce, NADHmean, NADmean, LDHmean)

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
            ncells[1,i] = 4;
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
        "ScFms" => ScFms
    
    );

    return data_multi

end

function restructDatInfMM_RealExpDat(paths, conce, LDHmean)

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

    Y0us = Array{Float64,2}(undef, 9, m);
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
        ScFms[i] = DataAll[i][1,2]./((conce*1000)*0.12)[1];
        # Pout, Pin, LDH, LDHPyr, Lac, LacHyper, PyrHyper, Obs_LacHyper, Obs_PyrHyper = u;
        Y0us[:,i] = [DataAll[i][1,2]/(ScFms[i]*0.12), 0,LDHmean, 0, 0, 0, DataAll[i][1,2]/(ScFms[i]*0.12), 0, DataAll[i][1,2]]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
            ncells[1,i] = 4;
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
        "ScFms" => ScFms
    
    );

    return data_multi

end

function restructDatInfMM(paths, LDHmean)

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

    Y0us = Array{Float64,2}(undef, 5, m);
    tml = maximum([length(DataAll[i][1,1]:1:DataAll[i][end,1]) for i in 1:m]);
    ts = Array{Float64,2}(undef, tml, m);
    itp = Array{Float64,1}(undef, m);

    ms = collect(2:2:obser*2);


    for i in 1:m
        sts[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        sts2[1:length(DataAll[i][:,1]),i] = convert.(Int, DataAll[i][1:end,1]);
        ts[1:length(DataAll[i][1,1]:1:DataAll[i][end,1]),i] = convert.(Int, DataAll[i][1,1]:1:DataAll[i][end,1]);
        ts[1,i] = 1e-20;
        for j in 1:obser
            Means[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]];
            Erros[1:length(DataAll[i][:,1]), i, j] = DataAll[i][1:end,ms[j]+1];
        end
        # Pout, Pin, LDH, LDHPyr, Lac, LacHyper, PyrHyper, Obs_LacHyper, Obs_PyrHyper = u;
        Y0us[:,i] = [DataAll[i][1,2], 0,LDHmean, 0, 0]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
            ncells[1,i] = 4;
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
    );

    return data_multi

end


function restructDatInfCompRep_RealExpDat_FixScF(paths, conce, NADHmean, NADmean, LDHmean, ScFEx)

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
        ScFms[i] = ScFEx;
        Y0us[:,i] = [DataAll[i][1,2]/(ScFms[i]*0.12), 0,0, NADHmean, NADmean, LDHmean, 0, 0, 0, 0, 0, DataAll[i][1,2]/(ScFms[i]*0.12), 0, DataAll[i][1,2]]; # Will need a better way for this!!! Perhaps introduce Y0 as an input to the function

        itp[i] = 0;
    end

    ncells = zeros(1, m);
    for i in 1:m
            ncells[1,i] = 4;
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
        "ScFms" => ScFms
    
    );

    return data_multi

end