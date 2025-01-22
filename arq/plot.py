import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Definir o diretório de trabalho (substitua pelo diretório desejado)
import os
#os.chdir('caminho/do/diretório')

# Criar a lista de arquivos e ângulos correspondentes
files = []
angles = []
for i in range(0, 10):  # Ajustado para corresponder ao loop original
    file_name = f"leopardv30{str(i).zfill(2)}.dat"
    if os.path.isfile(file_name):
        angle = (i % 10) * 10
        if i == 9 or i == 19 or i == 29:
            angle = 88
        print(f"Arquivo encontrado: {file_name} (Ângulo: {angle}°)")
        files.append(file_name)
        angles.append(angle)
    else:
        print(f"Aviso: O arquivo {file_name} não existe.")

# Inicializar listas para armazenar os dados
data = []

# Carregar os dados dos arquivos .dat
for file in files:
    try:
        d = np.loadtxt(file)
        if d.size > 0 and d.shape[1] >= 3:
            data.append(d)
        else:
            print(f"Aviso: O arquivo {file} está vazio ou não contém dados suficientes.")
    except Exception as e:
        print(f"Erro ao ler o arquivo {file}: {e}")

# Verifique se há dados para plotar
if len(data) == 0:
    raise ValueError("Nenhum dado válido foi carregado.")

# Inicializar o gráfico com limites dos eixos definidos
x_limits = (0, 5.0)
y_limits = (-1.0, 10.0)

# Cores diferentes para cada faixa de arquivos
colors = cm.get_cmap('tab10', len(data))

# Inicializar conjunto para rastrear rótulos de legenda já adicionados
legend_labels = set()

# Configurar a plotagem
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Extrair as colunas e plotar os dados
for i, (d, angle) in enumerate(zip(data, angles)):
    x = d[:, 0]  # Primeira coluna para o eixo X
    y_real = d[:, 1]  # Segunda coluna - Parte real
    y_imag = d[:, 2]  # Terceira coluna - Parte imaginária
    color = colors(i)  # Escolhe a cor com base na faixa de arquivos
    
    # Adicionar rótulos à legenda uma vez por ângulo
    if angle in legend_labels:
        ax1.plot(x, y_real, linestyle='-', linewidth=2, color=color)
        ax2.plot(x, y_imag, linestyle='--', linewidth=2, color=color)
    else:
        ax1.plot(x, y_real, linestyle='-', linewidth=2, color=color, label=f'Real {angle}°')
        ax2.plot(x, y_imag, linestyle='--', linewidth=2, color=color, label=f'Imag {angle}°')
        legend_labels.add(angle)

# Configurar os rótulos dos eixos e limites
ax1.set_xlabel('kc/ω_{pi}')
ax1.set_ylabel('ω_r/Ω_i')
ax1.set_xlim(x_limits)
ax1.set_ylim(y_limits)
ax1.legend(fontsize=10, loc='upper right')
ax1.set_title('Gráfico da Parte Real')

ax2.set_xlabel('kc/ω_{pi}')
ax2.set_ylabel('ω_{im}/Ω_i')
ax2.set_xlim(x_limits)
ax2.set_ylim(y_limits)
ax2.legend(fontsize=10, loc='upper right')
ax2.set_title('Gráfico da Parte Imaginária')

# Ajustar layout
plt.tight_layout()

# Salvar o gráfico em um arquivo
plt.savefig("n=2_it=100_re_im.png")
plt.show()
