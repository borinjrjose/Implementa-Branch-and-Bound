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

function colorir_heuristica(nro_nos, qnt_cores, cor_no)
    melhor_solucao = Array{Int, 2}(undef, nro_nos, qnt_cores)
    quant_no_cor = zeros(Int, qnt_cores)
    cor_mais_abundante = -1

    for cor in cor_no
        quant_no_cor[cor] = quant_no_cor[cor] + 1
        if cor_mais_abundante == -1 || quant_no_cor[cor] > quant_no_cor[cor_mais_abundante]
            cor_mais_abundante = cor
        end
    end

    for i in 1:nro_nos
        for j in 1:qnt_cores
            if j == cor_mais_abundante
                melhor_solucao[i, j] = 1
            else
                melhor_solucao[i, j] = 0
            end
        end
    end

    limitante_superior = 0
    for cor in 1:qnt_cores
        if cor != cor_mais_abundante
            limitante_superior = limitante_superior + quant_no_cor[cor]
        end
    end

    return melhor_solucao, limitante_superior
end

function criar_modelo_dicionario_base(nro_nos, qnt_cores, cor_no)
    modelo = Model(Gurobi.Optimizer)

    x = Dict()
    for i in 1:nro_nos
        for j in 1:qnt_cores
            x[i,j] = @variable(modelo, lower_bound=0)
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
                    add_to_expression!(restricao_cor_conexa, 1, x[q,j])
                    add_to_expression!(restricao_cor_conexa, -1, x[k,j])

                    @constraint(modelo, restricao_cor_conexa <= 1)
                end
            end
        end
    end

    return modelo, x
end

function encontrar_variavel_branch(x, nro_nos, qnt_cores)
    for i in 1:nro_nos
        for j in 1:qnt_cores
            if abs(value(x[i, j])-round(value(x[i, j]))) > ϵ
                return x[i, j]
            end
        end
    end
end

function nova_melhor_solucao(melhor_solucao, x, nro_nos, qnt_cores)
    for i in 1:nro_nos
        for j in 1:qnt_cores
            melhor_solucao[i, j] = value(x[i, j])
        end
    end
end

function resolve_sem_gurobi(nro_nos, qnt_cores, cor_no)
    melhor_solucao, limitante_superior = colorir_heuristica(nro_nos, qnt_cores, cor_no)

    modelo, x = criar_modelo_dicionario_base(nro_nos, qnt_cores, cor_no)

    optimize!(modelo)

    if termination_status(modelo) == MathOptInterface.OPTIMAL
        limitante_inferior = getobjectivevalue(modelo)

        x_branch = encontrar_variavel_branch(x, nro_nos, qnt_cores)

        if x_branch == nothing
            limitante_superior = limitante_inferior
            nova_melhor_solucao(melhor_solucao, x, nro_nos, qnt_cores)
        else
            cons = @constraint(modelo, x_branch == 0)
            limitante_superior = branch_and_bound(modelo, x, nro_nos, qnt_cores, limitante_superior, melhor_solucao)
            delete(modelo, cons)

            cons = @constraint(modelo, x_branch == 1)
            limitante_superior = branch_and_bound(modelo, x, nro_nos, qnt_cores, limitante_superior, melhor_solucao)
            delete(modelo, cons)
        end

        println(melhor_solucao)
    else
        println(string("Erro: solver não encontrou solução ótima. Status = ", termination_status(modelo)))
    end
end

function branch_and_bound(modelo, x, nro_nos, qnt_cores, limitante_superior, melhor_solucao)
    optimize!(modelo)

    if termination_status(modelo) == MathOptInterface.OPTIMAL
        x_branch = encontrar_variavel_branch(x, nro_nos, qnt_cores)

        limitante_inferior_atual = getobjectivevalue(modelo)

        if limitante_inferior_atual >= (limitante_superior-ϵ)
            return limitante_superior
        else
            x_branch = encontrar_variavel_branch(x, nro_nos, qnt_cores)

            if x_branch == nothing
                limitante_superior = limitante_inferior_atual
                nova_melhor_solucao(melhor_solucao, x, nro_nos, qnt_cores)
                return limitante_superior
            else
                cons = @constraint(modelo, x_branch == 0)
                limitante_superior = branch_and_bound(modelo, x, nro_nos, qnt_cores, limitante_superior, melhor_solucao)
                delete(modelo, cons)

                cons = @constraint(modelo, x_branch == 1)
                limitante_superior = branch_and_bound(modelo, x, nro_nos, qnt_cores, limitante_superior, melhor_solucao)
                delete(modelo, cons)

                return limitante_superior
            end
        end
    else
        return limitante_superior
    end
end

function main()
    instancia_problema = "instancias/rand_10_2.txt"
    entrada = processar_entrada(instancia_problema)

    resolve_sem_gurobi(entrada[1], entrada[2], entrada[3])
end

main()
