using DelimitedFiles
using Plots

# Verifique o diretório de trabalho atual
println("Diretório de Trabalho Atual: ", pwd())

# Criar a lista de arquivos e ângulos correspondentes
files = []
angles = []
nroots = 3
for i in 0:29
    file_name = "leopardv30$(lpad(i, 2, '0')).dat"
    if isfile(file_name)
        angle = (i % 10) * 10  # Cálculo do ângulo para arquivos 00 a 09
        if i == 9 || i == 19 || i == 29
            angle = 88  # Ajuste para o ângulo 90º
        end
        println("Arquivo encontrado: $file_name (Ângulo: $angle\u00B0)")
        push!(files, file_name)
        push!(angles, angle)
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
        if size(d, 1) > 0 && size(d, 2) >= 2*nroots+1  # Verificando se a oitava coluna existe
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
x_limits = (0, 3.0)   # Defina os limites do eixo X
y_limits = (0, 2.5)  # Defina os limites do eixo Y

plot(title="Gráfico de Dados Múltiplos", xlabel="Eixo X", ylabel="Eixo Y", xlims=x_limits, ylims=y_limits)

# Estilos de linha diferentes suportados pela backend GR
line_styles = [:solid, :dash, :dashdot, :dashdotdot, :dot]

# Cores diferentes para cada faixa de arquivos
colors = ["green", "red", "blue"]

# Inicializar conjunto para rastrear rótulos de legenda já adicionados
legend_labels = Set{Int}()

# Extrair as colunas e plotar os dados
nreal=nroots-1
ncomplex=nroots
for (i, (d, angle)) in enumerate(zip(data, angles))
    x = d[:, 1]  # Primeira coluna (indexação começa em 1)
    y1 = d[:, nreal]  # raízes reais
    y2 = d[:, ncomplex]  # raízes complexas
    linestyle = line_styles[mod((angle ÷ 10), length(line_styles)) + 1]  # Escolhe o estilo da linha
    color = colors[div(i - 1, 10) + 1]  # Escolhe a cor com base na faixa de arquivos
    # Adicionar rótulo apenas se não estiver no conjunto de rótulos
    if angle in legend_labels
        plot!(x, y1, linestyle=linestyle, linewidth=2, color=color, label="")
    else
        plot!(x, y1, label="$angle\u00B0", linestyle=linestyle, linewidth=2, color=color)
        push!(legend_labels, angle)
    end
end

# Configurar a legenda para ser preta
plot!(legendfontsize=10, legend_foreground_color=:black, legend_background_color=:white, legend=:outertopright, legendcolor=:black)

# Salvar o gráfico em um arquivo
savefig("n=$nroots it=100.png")
display(plot)
