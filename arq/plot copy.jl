using DelimitedFiles
using Plots

# Verifique o diretório de trabalho atual
println("Diretório de Trabalho Atual: ", pwd())

# Criar a lista de arquivos
files = []
for i in 1:30
    file_name = "leopardv30$(lpad(i, 2, '0')).dat"
    if isfile(file_name)
        println("Arquivo encontrado: $file_name")
        push!(files, file_name)
    else
        println("Aviso: O arquivo $file_name não existe.")
    end
end

# Inicializar arrays para armazenar os dados
data = []

# Carregar os dados dos arquivos .dat
for file in files
    try
        d = readdlm(file)
        if size(d, 1) > 0 && size(d, 2) >= 3
            push!(data, d)
        else
            println("Aviso: O arquivo $file está vazio ou não contém dados suficientes.")
        end
    catch e
        println("Erro ao ler o arquivo $file: ", e)
    end
end

# Verifique se há dados para plotar
if length(data) == 0
    error("Nenhum dado válido foi carregado.")
end

# Inicializar o gráfico
x_limits = (0, 0.6)   # Defina os limites do eixo X
y_limits = (-0.05, 0.1)  # Defina os limites do eixo Y

plot(title="Gráfico de Dados Múltiplos", xlabel="Eixo X", ylabel="Eixo Y", xlims=x_limits, ylims=y_limits)

# Extrair as colunas e plotar os dados
for (i, d) in enumerate(data)
    x = d[:, 1]  # Primeira coluna (indexação começa em 1)
    y = d[:, 3]  # Terceira coluna (indexação começa em 1)
    plot!(x, y, label="Dados $i", marker=:circle)
end

# Salvar o gráfico em um arquivo
savefig("grafico_multiplos_dados.png")
display(plot)
