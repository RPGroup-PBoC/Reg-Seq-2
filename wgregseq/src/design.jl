using DataFrames, CSV, wgregseq, HTTP, LightXML

struct gene
    name::String
    ID::String
end

function get_ecocyc(genes::String)
    return get_ecocyc(String[genes])

end


"""
    get_ecocyc(genes::Array{String, 1})

Get information about gene from Ecocyc.
"""
function get_ecocyc(genes::Array{String, 1})
    path = dirname(dirname(pathof(wgregseq))) * "/data/genes_MG1655.txt"
    df = CSV.read(path, DataFrame)
    mask = [gene in genes for gene in df[!, "Gene Name"]]
    sub_df = df[mask, :]
    gene_list = gene[]
    xml_list = []
    for i in 1:nrow(sub_df)
        name, ID = sub_df[i, ["Gene Name", "Object ID"]]
        r = HTTP.get("https://websvc.biocyc.org/apixml?fn=transcription-units-of-gene&id=ECOLI:$ID&detail=full")
        xdoc = parse_string(String(copy(r.body)))
        
        push!(xml_list, xdoc)
        push!(gene_list, gene(name, ID))
        sleep(1)
    end
    return xml_list    
end