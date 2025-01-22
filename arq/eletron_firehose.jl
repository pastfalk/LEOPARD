#import Pkg
#Pkg.add("LaTeXStrings")
using DelimitedFiles
using Plots
using LaTeXStrings  # Para usar LaTeX nos gráficos

# Configurar o backend GR
gr()

# Carregar o arquivo
data = readdlm("leopardv3000.dat")

# Extrair colunas relevantes
x = data[:, 1]           # Primeira coluna (kc/ωpi)
real_part = data[:, 6]   # Parte real (ωr/Ωi)
imag_part = data[:, 7]   # Parte imaginária (γ/Ωi)

# Configurar layout para gráficos sobrepostos
l = @layout [a; b]

# Gráfico da parte real
p1 = plot(x, real_part, label="", xlabel="", ylabel=L"ω_r/Ω_i", linecolor=:black, linewidth=2, legend=false)
plot!(p1, x, real_part, linestyle=:dash, linecolor=:black)

# Gráfico da parte imaginária
p2 = plot(x, imag_part, label="", xlabel=L"kc/ω_{pi}", ylabel=L"ω_{im}/Ω_i", linecolor=:black, linewidth=2, legend=false)
plot!(p2, x, imag_part, linestyle=:dash, linecolor=:black)

# Ajustar limites dos eixos
x_limits = (0, 5.0)   # Limites do eixo X
y1_limits = (0, 10.0) # Limites do eixo Y para o gráfico superior
y2_limits = (-0.5, 0.5)  # Limites do eixo Y para o gráfico inferior

plot!(p1, xlims=x_limits, ylims=y1_limits)
plot!(p2, xlims=x_limits, ylims=y2_limits)

# Combinar os gráficos
final_plot = plot(p1, p2, layout=l, size=(800, 600))

# Exibir o gráfico
gui()

# Salvar o gráfico com resolução ajustada
savefig(final_plot, "electron_firehose_instability.png")