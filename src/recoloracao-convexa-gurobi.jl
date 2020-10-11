using JuMP
using Gurobi
using MathOptInterface

ϵ = 0.000001

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

function resolve_com_gurobi(nro_nos, qnt_cores, cor_no)
    modelo = Model(Gurobi.Optimizer)

    x = Dict()
    for i in 1:nro_nos
        for j in 1:qnt_cores
            x[i,j] = @variable(modelo, binary=true)
        end
    end

    qnt_troca_cores = AffExpr()
    for i in 1:nro_nos
        for j in 1:qnt_cores
            if cor_no[i] != j
                add_to_expression!(qnt_troca_cores, 1, x[i,j])
            end
        end
    end

    @objective(modelo, Min, qnt_troca_cores)

    for i in 1:nro_nos
        max_cor_no = AffExpr()
        for j in 1:qnt_cores
            add_to_expression!(max_cor_no, 1, x[i,j])
        end

        @constraint(modelo, max_cor_no == 1)
    end

    for j in 1:qnt_cores
        for q in 3:nro_nos
            for p in 1:(q-2)
                for k in (p+1):(q-1)
                    restricao_cor_conexa = AffExpr()
                    add_to_expression!(restricao_cor_conexa, 1, x[p,j])
                    add_to_expression!(restricao_cor_conexa, -1, x[k,j])
                    add_to_expression!(restricao_cor_conexa, 1, x[q,j])

                    @constraint(modelo, restricao_cor_conexa <= 1)
                end
            end
        end
    end

    optimize!(modelo)

    if termination_status(modelo) == MathOptInterface.OPTIMAL
        imprime_cor_grafo(x, nro_nos, qnt_cores)
    else
        println(string("Erro: solver não encontrou solução ótima. Status = ", termination_status(modelo)))
    end
end

function imprime_cor_grafo(x, nro_nos, qnt_cores)
    println()
    for i in 1:nro_nos
        for j in 1:qnt_cores
            if value(x[i, j]) > (1-ϵ)
                print(string(j, " "))
            end
        end
    end
end

function main()
    instancia_problema = "instancias/rand_10_2.txt"
    entrada = processar_entrada(instancia_problema)

    resolve_com_gurobi(entrada[1], entrada[2], entrada[3])
end

main()
