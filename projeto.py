# O Projeto: um panorama da família Flaviviridae
#
# Leia o enunciado completo no README (seção "O Projeto")
#
# A ideia é construir UMA tabela (pandas) descrevendo os vírus e, a partir dela,
# tirar duas conclusões:
#   - o conteúdo GC é aleatório? (Parte 2)
#   - quão grande é a proteína de cada vírus? (Parte 3)
#
# Vá preenchendo as partes abaixo, uma de cada vez.
# Obs: Se preferir fazer esse processo num jupyter notebook, sem problemas!! Fica a critério do grupo

import pandas as pd

pd.set_option('display.max_columns', None)   # mostra todas as colunas
pd.set_option('display.width', None)         # não quebra linha por largura do terminal
pd.set_option('display.max_colwidth', 30)

from bio.ler_fasta import ler_fasta
from bio.sequencia import (
    traduzir,
    calcular_percentual_gc,
    encontrar_inicio,
)


# ------------------------------------------------------------------
# Parte 1 — Monte a tabela
# ------------------------------------------------------------------
organismos = ler_fasta("arquivos/Flaviviridae-genomes.fasta")
df = pd.DataFrame(organismos)
#print(df.head())




# ------------------------------------------------------------------
# Parte 2 — O conteúdo GC é aleatório?
# ------------------------------------------------------------------
# 1) crie a coluna "gc" com df["sequencia"].apply(calcular_percentual_gc)
df["gc"] = df["sequencia"].apply(calcular_percentual_gc)
#print(df.head())

# 2) mostre os 10 maiores e os 10 menores GC (com o nome!) -> usar função sort_values do pandas
df_ordenado = df.sort_values(by='gc')
dez_menores_gc = df_ordenado[["nome", "gc"]].head(10)
dez_maiores_gc = df_ordenado[["nome", "gc"]].tail(10)
print(f'Os 10 vírus com menor conteúdo GC em suas sequência são: {dez_menores_gc}')
print(f'Os 10 vírus com maior conteúdo GC em suas sequência são: {dez_maiores_gc}')

# 3) escreva sua conclusão sobre o padrão que observou
'''É possível notar que os vírus com menor conteúdo GC são pertencentes ao gênero Pestivirus.
Enquanto isso, os vírus com maior conteúdo GC são pertencentes ao gênero Pegivirus e Hepacivirus.'''

# Parte 3 — Encontre a proteína (a poliproteína viral)
# ------------------------------------------------------------------
# 1) coluna "proteina": traduzir(encontrar_inicio(seq), parar=True)
df["proteina"] = df["sequencia"].apply(
    lambda sequencia: traduzir(encontrar_inicio(sequencia), parar=True)
)
# 2) coluna "tamanho_proteina": len da proteína
df["tamanho_proteina"] = df["proteina"].apply(len)
# 3) coluna "cobertura": (tamanho_proteina * 3) / tamanho
df["tamanho"] = df["sequencia"].apply(len)
df["cobertura"] = (df["tamanho_proteina"] * 3) / df["tamanho"]
print(df.head(10))
# 4) escreva sua conclusão (qual a cobertura típica? faz sentido ser 1 poliproteína?)
'''
A cobertura encontrada apresentou valores predominantemente altos, 
com mediana de aproximadamente 0,93 e 75% dos vírus apresentando cobertura acima de 0,95.
Isso indica que a região traduzida corresponde à maior parte do genoma viral.
Esse resultado é esperado para vírus da família Flaviviridae, 
que tem uma única sequência codificante longa responsável pela produção de uma poliproteína, 
posteriormente processada em diferentes proteínas virais.
'''


# ------------------------------------------------------------------
# Parte 4 — Salve o resultado
# ------------------------------------------------------------------
# 1) filtre os vírus com gc > 0.5 (quantos são?)
df_gc_alto = df[df["gc"] > 0.5]
print(f"Vírus com GC > 0.5: {len(df_gc_alto)} de {len(df)}")
print(df_gc_alto[["nome", "gc"]])

# 2) df.to_csv("resultado.csv", index=False)
df.to_csv("resultado.csv", index=False)
