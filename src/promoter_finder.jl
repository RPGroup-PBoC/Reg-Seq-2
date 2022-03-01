using BioSequences, GLM, DataFrames, CSV

####################################
## Following Code is adapted from ##
## La Fleur et al. 2022           ##
####################################
# k and BETA for La Fleur dataset
LOGK   = -2.80271176
BETA    = 0.81632623

# Dictionaries
groove_access = Dict(
    "CG" => 43,
    "CA"=>  42, "TG" => 42, #revcomp
    "GG"=>  42, "CC" => 42,  #revcomp
    "GC"=>  25, 
    "GA"=>  22, "TC" => 22, #revcomp
    "TA"=>  14,
    "AG"=>  9, "CT" => 9, #revcomp
    "AA"=>  5, "TT" => 5, #revcomp
    "AC"=>  4, "GT" => 4, #revcomp
    "AT"=>  0
)

stacking_dict = Dict(
    "AA" =>-2.41,      
    "AC"=>-2.06, "CA" => -2.06, #reverse   
    "AG"=>-2.54, "GA" => -2.54, #reverse   
    "AT"=>-2.25, "TA" => -2.25, #reverse    
    "CC"=>-1.93,    
    "CG"=>-2.24, "GC" => -2.24, #reverse   
    "CT"=>-2.01, "TC" => -2.01, #reverse  
    "GG"=>-2.66,   
    "GT"=>-2.46, 
    "TG" => -2.46, #reverse   
    "TT"=>-1.82
)

average_twist = Dict(
    "AA" => 35.6, "TT" => 35.6, #revcomp
    "AT"=>  32.1, "AT" => 32.1, #revcomp
    "TA"=>  35.3, "TA" => 35.3, #revcomp
    "GG"=>  33.65, "CC" => 33.65, #revcomp
    "GC"=>  40.2,
    "CG"=>  29.9,
    "AG"=>  27.9, "CT" => 27.9, #revcomp
    "GA"=>  36.8, "TC" => 36.8, #revcomp
    "AC"=>  34.0, "GT" => 34.0, #revcomp
    "CA"=>  34.4, "TG" => 34.4 #revcomp
)

# Sequence dependence of DNA bending rigidity, Geggier and Vologodskii
persistence = Dict(
    "AA" => 50.4, "TT" => 50.4, #revcomp
    "AC"=>  55.4, "GT" => 55.4, #revcomp
    "AG"=>  51.0, "CT" => 51.0, #revcomp
    "AT"=>  40.9, "AT" => 40.9, #revcomp 
    "CA"=>  46.7, "TG" => 46.7, #revcomp
    "CC"=>  41.7, "GG" => 41.7, #revcomp
    "CG"=>  56.0,
    "GA"=>  54.4, "TC" => 54.4, #revcomp
    "GC"=>  44.6,
    "TA"=>  44.7,
)

RNA_DNA_hybrids = Dict(
    "TT" =>  -1.0, "AA" => -1.0,  #revcomp
    "TG"=>   -2.1, "CA" => -2.1,  #revcomp
    "TC"=>   -1.8, "GA" => -1.8,  #revcomp
    "TA"=>   -0.9,
    "GT"=>   -0.9, "AC" => -0.9,  #revcomp
    "GG"=>   -2.1, "CC" => -2.1,  #revcomp
    "GC"=>   -1.7,
    "GA"=>   -0.9, "TC" => -0.9,  #revcomp
    "CT"=>   -1.3, "AG" => -1.3,  #revcomp
    "CG"=>   -2.7,
    "CC"=>   -2.9, "GG" => -2.9,  #revcomp
    "CA"=>   -1.1, "TG" => -1.1,  #revcomp
    "AT"=>   -0.6,
    "AG"=>   -1.5, "CT" => -1.5,  #revcomp
    "AC"=>   -1.6, "GT" => -1.6,  #revcomp
    "AA"=>   -0.2, "TT" => -0.2   #revcomp
)

DNA_DNA_hybrids = Dict( 
    "AA"=>   -1.00, "TT" => -1.00, #revcomp
    "AT"=>   -0.88,
    "TA"=>   -0.58,
    "CA"=>   -1.45, "TG" => -1.45, #revcomp
    "GT"=>   -1.44, "AC" => -1.44, #revcomp
    "CT"=>   -1.28, "AG" => -1.28, #revcomp
    "GA"=>   -1.30, "TC" => -1.30, #revcomp
    "CG"=>   -2.17,
    "GC"=>   -2.24,
    "GG"=>   -1.42, "CC" => -1.42  #revcomp
)


function get_model_params()
    dir = @__DIR__
    home_dir = joinpath(split(dir, '/')[1:end-1]...)
    # Initialize model and matrices
    layer1 = CSV.read("/$home_dir/data/salis/free_energy_coeffs.txt", DataFrame, header=false)[!, "Column1"]

    inters = open("/$home_dir/data/salis/model_intercept.txt") do file
        read(file, String)
    end
    inters = parse(Float64, inters)
    return layer1, inters

end


all_kmers(k) = vcat(join.(collect(Iterators.product(fill([DNA_A, DNA_C, DNA_G, DNA_T], k)...)))...) |> sort


### This function calcs the hamming dist between 2 sequences of = len
function hamming(seq1, seq2)
    if length(seq1) != length(seq2)
        throw(ArgumentError("Sequences are not of equal length! Sequence 1: $(length(seq1)), Sequence 2: $(length(seq2))"))
    end
    distance = sum([x1 != x2 for (x1, x2) in zip(collect(seq1), collect(seq2))])
    return distance
end


function get_rsquared(TXpred, TXempirical)
    ols = lm(@formula(A  ~ B), DataFrame(A=TXpred, B=TXempirical))
    return r2(ols)
end


function AT_content(seq)
    ATcont = count(x -> x == DNA_A, collect(seq)) + count(x -> x == DNA_T, collect(seq))
   
    return ATcont / length(seq) * 100
end


# Functions to calculate DNA or RNA properties
function calc_stacking_and_twist_energy(seq)
    #find stacking free enegy of spacer
    spacer_sfe  = 0
    spacer_twist = 0
    for j in range(1, stop=length(seq), step=2)
        mer = seq[j:j+1]
        spacer_sfe += stacking_dict[string(mer)]
        spacer_twist += average_twist[string(mer)]
    end
    return spacer_sfe, spacer_twist
end


# local DNA element rigidity
function calc_rigidity(seq)
    rigidity = 0
    for m in range(1, stop=length(seq)-1, step=2)
        mer = seq[m:m+1]
        rigidity += persistence[string(mer)]
    end
    rigidity = rigidity / length(seq)
    return rigidity
end


function calc_groove_width(seq)
    groove_width = 0
    for l in range(1, stop=length(seq)-1, step=2)
        mer1 = seq[l:l+1]
        groove_width += groove_access[string(mer1)]
    end
    return groove_width
end


function calc_DNA_RNA_hybrid_energy(seq)
    dg_dna = 0
    dg_rna = 0
    for j in range(1, stop=length(seq[1:15])-1, step=2)
         # 15 is the length of long abortive transcripts
        mer = seq[j:j+1]
        dg_dna += DNA_DNA_hybrids[string(mer)]
        dg_rna += RNA_DNA_hybrids[string(mer)]
    end
    dg_hybrid = dg_dna - dg_rna
    return dg_dna, dg_rna, dg_hybrid
end


function get_matrices(coeffs)
    three_mer = LongDNA{4}.(all_kmers(3))
    two_mer = LongDNA{4}.(all_kmers(2))
    #Extract dG values from model coefficients
    dg10_0 = Dict(three_mer .=> coeffs[1:64])
    dg10_3 = Dict(three_mer .=> coeffs[65:128])
    dg35_0 = Dict(three_mer .=> coeffs[129:192])
    dg35_3 = Dict(three_mer .=> coeffs[193:256])
    dmers = Dict(three_mer .=> coeffs[257:256+64])
    x10mers = Dict(two_mer .=> coeffs[256+64+1:256+64+16])
    spacers = Dict(string.(collect(16:18)) .=> [256+64+16+1:256+64+16+3])

    return dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers
end


function linear_free_energy_model(
    UP, 
    h35, 
    spacer, 
    h10, 
    disc, 
    ITR, 
    dg10_0, 
    dg10_3, 
    dg35_0, 
    dg35_3, 
    dmers, 
    x10mers,
    coeffs, 
    inters
    )

    prox_UP = UP[end-Int(floor(length(UP)/2))+1:end]
    dist_UP = UP[1:Int(floor(length(UP)/2))]

    # CATEGORICAL FEATURES
    ext10           = spacer[end-2:end-1] # TGN motif, contacts sigma
    hex10_0         = h10[1:3]
    hex10_3         = h10[4:end]
    hex35_0         = h35[1:3]
    hex35_3         = h35[4:end]
    disc_first_3mer = disc[1:3]
    spacer_length   = string(length(spacer))

    # NUMERICAL FEATURES
    dg_dna, dg_rna, dg_ITR     = calc_DNA_RNA_hybrid_energy(ITR) # calc R-loop strength
    rigidity                   = calc_rigidity(join([UP, h35, spacer[1:14]]))

    width_proxy_prox = calc_groove_width(prox_UP)
    width_proxy_dist = calc_groove_width(dist_UP)

    # NORMALIZE NUMERICAL FEATURES BY MAX IN TRAINING SET
    numericals         = [width_proxy_dist, width_proxy_prox, dg_ITR, rigidity]
    normalizing_values = [256.0, 255.0, 4.300000000000002, 25.780434782608694]
    numerical_coefs    = coeffs[end-3:end]
    normald            = numericals ./ normalizing_values
    dg_numerical       = normald .* numerical_coefs

    dg10      = dg10_0[hex10_0] + dg10_3[hex10_3]
    dg35      = dg35_0[hex35_0] + dg35_3[hex35_3]
    dg_disc   = dmers[disc_first_3mer]
    dg_ITR    = dg_numerical[end-1]
    dg_ext10  = x10mers[ext10]
    x = parse(Float64, spacer_length)
    dg_spacer = 0.1463*x^2 - 4.9113*x + 41.119

    dg_UP        = dg_numerical[1] + dg_numerical[2] + dg_numerical[end]
    dG_apparent  = (dg10 + dg35 + dg_disc + dg_ITR + dg_ext10 + dg_spacer + dg_UP + inters[1] - LOGK)/BETA
    dG_total     = dg10 + dg35 + dg_disc + dg_ITR + dg_ext10 + dg_spacer + dg_UP + inters[1]

    return dG_total, dG_apparent, dg10, dg35, dg_disc, dg_ITR, dg_ext10, dg_spacer, dg_UP
end




struct Promoter_Calculator
    model
    inters
    dg10_0
    dg10_3
    dg35_0
    dg35_3
    dmers
    x10mers
    spacers
    K
    BETA
    function Promoter_Calculator()
        layer1, inters = get_model_params()
        dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers = get_matrices(layer1)
        K = 42.00000
        BETA = 1.636217004872062
        return new(layer1, inters, dg10_0, dg10_3, dg35_0, dg35_3, dmers, x10mers, spacers, K, BETA)
    end
end



function (p::Promoter_Calculator)(sequence::BioSequences.LongDNA, TSS_range::Tuple{Int64, Int64}=(-1, -1))
   
    if TSS_range == (-1, -1) 
        TSS_range = (1, length(sequence))
    end

    TSS_range_rev = (length(sequence) - TSS_range[2]+1, length(sequence) - TSS_range[1]+1)
    fwd_sequence = sequence
    rev_sequence = reverse_complement(sequence)
    Forward_Min_States, Forward_All_States = predict(p, fwd_sequence, TSS_range)
    Reverse_Min_States_Temp, Reverse_All_States_Temp = predict(p, rev_sequence, TSS_range_rev)

    Reverse_Min_States = Dict()
    Reverse_All_States = Dict()

    for TSS in Reverse_Min_States_Temp |> keys |> collect
        Reverse_Min_States[length(sequence) - TSS] = Reverse_Min_States_Temp[TSS]
        Reverse_All_States[length(sequence) - TSS] = Reverse_All_States_Temp[TSS]
    end
    Forward_Predictions_per_TSS = Forward_Min_States
    Reverse_Predictions_per_TSS = Reverse_Min_States
    output = Dict(
        "K" => p.K,
        "beta" => p.BETA,
        "sequence" => sequence,
        "TSS_range" => TSS_range,
        "Forward_Predictions_per_TSS" => Forward_Predictions_per_TSS,
        "Reverse_Predictions_per_TSS" => Reverse_Predictions_per_TSS
    )
    return output
end


function predict(p::Promoter_Calculator, sequence::BioSequences.LongDNA, TSS_range::Tuple{Int64, Int64})
    UPS_length = 24
    HEX35_length = 6
    UPS_HEX35_SPACER = 1
    SPACER_length_range = [15, 20]
    HEX10_length = 6
    DISC_length_range = [6, 10]
    ITR_length = 20

    MinPromoterSize = UPS_length + UPS_HEX35_SPACER + HEX35_length + SPACER_length_range[1] + HEX10_length + DISC_length_range[1] + ITR_length
    MaxPromoterSize = UPS_length + UPS_HEX35_SPACER + HEX35_length + SPACER_length_range[2] + HEX10_length + DISC_length_range[2] + ITR_length
    MinimumTSS = UPS_length + UPS_HEX35_SPACER + HEX35_length + SPACER_length_range[1] + HEX10_length + DISC_length_range[1]

    All_States = Dict()
    Min_States = Dict()
    for TSS in TSS_range[1]:TSS_range[2]
        All_States[TSS] = Dict()
        for DISC_length in DISC_length_range[1]:DISC_length_range[2]
            if (TSS - DISC_length >= 0) && (TSS + ITR_length <= length(sequence))
                tempdisc = sequence[TSS - DISC_length + 1: TSS]
                tempITR  = sequence[TSS + 1 : TSS + ITR_length]
                
                for SPACER_length in SPACER_length_range[1]:SPACER_length_range[2]
                    if TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_length - UPS_HEX35_SPACER >= 0
                        
                        temp10     = sequence[ TSS - DISC_length - HEX10_length + 1 : TSS - DISC_length ]
                        tempspacer = sequence[ TSS - DISC_length - HEX10_length - SPACER_length + 1 : TSS - DISC_length - HEX10_length]
                        temp35     = sequence[( TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length +1): TSS - DISC_length - HEX10_length - SPACER_length]
                        
                        tempUP     = sequence[ TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_length - UPS_HEX35_SPACER + 1:  TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_HEX35_SPACER]
                        
                        
                        dG_total, dG_apparent, dG_10, dG_35, dG_disc, dG_ITR, dG_ext10, dG_spacer, dG_UP = linear_free_energy_model(tempUP, temp35, tempspacer, temp10, tempdisc, tempITR, p.dg10_0, p.dg10_3, p.dg35_0, p.dg35_3, p.dmers, p.x10mers, p.model, p.inters)
                        
                        dG_bind  = dG_10 + dG_35 + dG_spacer + dG_ext10 + dG_UP

                        Tx_rate = p.K * exp(- p.BETA * dG_total)

                        result = Dict(
                            "promoter_sequence" => sequence[TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_length - UPS_HEX35_SPACER + 1: TSS + ITR_length],
                            "TSS" => TSS, 
                            "UP" => tempUP, 
                            "hex35" => temp35, 
                            "spacer" => tempspacer, 
                            "hex10" => temp10, 
                            "disc" => tempdisc, 
                            "ITR" => tempITR,
                            "dG_total" => dG_total, 
                            "dG_10" => dG_10, 
                            "dG_35" => dG_35, 
                            "dG_disc" => dG_disc, 
                            "dG_ITR" => dG_ITR, 
                            "dG_ext10"=> dG_ext10, 
                            "dG_spacer" => dG_spacer, 
                            "dG_UP" => dG_UP, 
                            "dG_bind" => dG_bind,
                            "Tx_rate" => Tx_rate,
                            "UP_position" => [TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_HEX35_SPACER - UPS_length + 1,
                                                       TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length - UPS_HEX35_SPACER],
                            "hex35_position" => [TSS - DISC_length - HEX10_length - SPACER_length - HEX35_length + 1, TSS - DISC_length - HEX10_length - SPACER_length],
                            "spacer_position" => [TSS - DISC_length - HEX10_length - SPACER_length + 1, TSS - DISC_length - HEX10_length],
                            "hex10_position" => [TSS - DISC_length - HEX10_length + 1, TSS - DISC_length],
                            "disc_position" => [TSS - DISC_length + 1, TSS]
                        )

                        All_States[TSS][(DISC_length, SPACER_length)] = result
                        if TSS in Min_States |> keys |> collect
                            if result["dG_total"] < Min_States[TSS]["dG_total"]  
                                Min_States[TSS] = result
                            end
                        else
                            Min_States[TSS] = result
                        end
                    end
                end
            end
        end
    end
    return (Min_States, All_States)
end

