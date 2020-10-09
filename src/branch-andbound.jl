using JuMP
using Gurobi
using MathOptInterface

function processar_entrada(nome_arq)
    linha_arq = []

    open(nome_arq) do arq
        linha_arq = readlines(arq)
    end

    nro_nos, qnt_cores = parse.(Int, split(linha_arq[1]))
    cor_no = []

    for i in 2:(nro_nos+1)
        push!(cor_no, parse(Int, linha_arq[i]))
    end

    return nro_nos, qnt_cores, cor_no
end

function main()
    instancia_problema = "instancias/rand_10_2.txt"
    entrada = processar_entrada(instancia_problema)
end

main()
