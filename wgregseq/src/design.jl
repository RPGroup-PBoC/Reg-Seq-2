using DataFrames, CSV, wgregseq, HTTP, LightXML

struct gene
    name::String
    TSS::Vector{Int64}
    direction::Vector{String}
end

function get_ecocyc(genes::String)
    return get_ecocyc(String[genes])
end


"""
    get_ecocyc(genes::Array{String, 1})

Get information about gene from Ecocyc.
"""
function get_ecocyc(genes::Vector{String})
    path = dirname(dirname(pathof(wgregseq))) * "/data/genes_MG1655.txt"
    df = CSV.read(path, DataFrame)
    mask = [gene in genes for gene in df[!, "Gene Name"]]
    sub_df = df[mask, :]
    gene_list = gene[]
    
    for i in 1:nrow(sub_df)
        name, ID = sub_df[i, ["Gene Name", "Object ID"]]
        r = HTTP.get("https://websvc.biocyc.org/apixml?fn=transcription-units-of-gene&id=ECOLI:$ID&detail=full")
        xdoc = parse_string(String(copy(r.body)))
        xroot = root(xdoc) 
        ces = xroot["Transcription-Unit"]
        
        tss_list = Int64[]
        direction_list = String[]
        
        for tu in ces
            TU_ID = attribute(tu, "ID")
            sleep(1.5)
            prom = HTTP.get("https://websvc.biocyc.org/apixml?fn=transcription-unit-promoter&id=$(TU_ID)&detail=full")
            prom_xdoc = parse_string(String(copy(prom.body)))
            prom_xroot = root(prom_xdoc)
            if length(prom_xroot["Promoter"]) > 0
                if length(prom_xroot["Promoter"][1]["absolute-plus-1-pos"]) > 0
                    push!(tss_list, parse(Int64, content(prom_xroot["Promoter"][1]["absolute-plus-1-pos"][1])))
                end
                if length(prom_xroot["Promoter"][1]["transcription-direction"]) > 0
                    push!(direction_list, content(prom_xroot["Promoter"][1]["transcription-direction"][1]))
                end
            end
            
        end
        
        push!(gene_list, gene(name, tss_list, direction_list))
        
    end
    return gene_list    
end