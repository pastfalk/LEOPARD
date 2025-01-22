import Pkg
Pkg.add("LaTeXStrings")
using DelimitedFiles
using Plots
using LaTeXStrings  # Para usar LaTeX nos gráficos

gr()
# Verifique o diretório de trabalho atual
#println("Diretório de Trabalho Atual: ", pwd())

# Criar a lista de arquivos e ângulos correspondentes
files = []
angles = []
for i in 0:0
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

# Inicializar o gráfico com limites dos eixos definidos
x_limits = (0, 5.0)   # Defina os limites do eixo X
y_limits = (-1., 10.0)  # Defina os limites do eixo Y
# Configurar layout para gráficos sobrepostos
l = @layout [a; b]

# Gráfico da parte real
p1 = plot(x, real_part, label="", xlabel="", ylabel=L"ω_r/Ω_i", linecolor=:black, linewidth=2, legend=false)
plot!(p1, x, real_part, linestyle=:dash, linecolor=:black)

# Gráfico da parte imaginária
p2 = plot(x, imag_part, label="", xlabel=L"kc/ω_{pi}", ylabel=L"ω_{im}/Ω_i", linecolor=:black, linewidth=2, legend=false)
plot!(p2, x, imag_part, linestyle=:dash, linecolor=:black)

plot(title="Gráfico de Dados Reais e Imaginários", xlabel="Eixo X", ylabel="Eixo Y", xlims=x_limits, ylims=y_limits)
# Cores diferentes para cada faixa de arquivos
colors = ["green", "red", "blue"]

# Inicializar conjunto para rastrear rótulos de legenda já adicionados
legend_labels = Set{Int}()

# Extrair as colunas e plotar os dados
for (i, (d, angle)) in enumerate(zip(data, angles))
    x = d[:, 1]  # Primeira coluna para o eixo X
    y_real = d[:, 2]  # Segunda coluna - Parte real
    y_imag = d[:, 3]  # Terceira coluna - Parte imaginária
    color = colors[div(i - 1, 10) + 1]  # Escolhe a cor com base na faixa de arquivos
    
    # Adicionar rótulos à legenda uma vez por ângulo
    if angle in legend_labels
        plot!(x, y_real, linestyle=:solid, linewidth=2, color=color, label="")
        plot!(x, y_imag, linestyle=:dash, linewidth=2, color=color, label="")
    else
        plot!(x, y_real, label="Real $angle\u00B0", linestyle=:solid, linewidth=2, color=color)
        plot!(x, y_imag, label="Imag $angle\u00B0", linestyle=:dash, linewidth=2, color=color)
        push!(legend_labels, angle)
    end
end

# Configurar a legenda
plot!(legendfontsize=10, legend_foreground_color=:black, legend_background_color=:white, legend=:outertopright)

# Salvar o gráfico em um arquivo
savefig("n=2 it=100_re_im.png")
display(plot)
