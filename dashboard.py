"""
dashboard.py
============

Funções para consultar o banco de variantes e montar o painel (dashboard)
com os principais indicadores.

Autora: Aline Machado
Licença: MIT
"""

from __future__ import annotations

import os
import re
import sqlite3

import matplotlib.pyplot as plt
import pandas as pd

CORES_CATEGORIA = {
    "Patogênica": "#C8553D",
    "Benigna": "#4C9F70",
    "Significado incerto": "#E0A93B",
    "Conflitante": "#7B5EA7",
    "Outra": "#9E9E9E",
}
ORDEM_CATEGORIAS = list(CORES_CATEGORIA)
CORES_GENE = {"BRCA1": "#3A6EA5", "BRCA2": "#C8553D"}


def carregar_consultas(caminho: str | None = None) -> dict[str, str]:
    """Lê o arquivo consultas.sql e devolve {nome_da_consulta: sql}."""
    if caminho is None:
        caminho = os.path.join(os.path.dirname(os.path.abspath(__file__)), "consultas.sql")
    with open(caminho, encoding="utf-8") as f:
        texto = f.read()
    blocos = re.split(r"^-- nome:\s*(\w+)\s*$", texto, flags=re.MULTILINE)
    return {nome: sql.strip().rstrip(";") for nome, sql in zip(blocos[1::2], blocos[2::2])}


CONSULTAS = carregar_consultas()


def consultar(conn: sqlite3.Connection, nome: str) -> pd.DataFrame:
    """Executa uma consulta nomeada de consultas.sql e devolve um DataFrame."""
    return pd.read_sql_query(CONSULTAS[nome], conn)


def _br(numero) -> str:
    """Formata números no padrão brasileiro: 1.234 e 12,5."""
    if isinstance(numero, float):
        return f"{numero:.1f}".replace(".", ",")
    return f"{numero:,}".replace(",", ".")


def _tabela_categorias(df: pd.DataFrame, indice: str) -> pd.DataFrame:
    tabela = df.pivot_table(index=indice, columns="categoria", values="variantes",
                            aggfunc="sum", fill_value=0)
    return tabela[[c for c in ORDEM_CATEGORIAS if c in tabela.columns]]


def painel(conn: sqlite3.Connection, salvar_em: str | None = None):
    """Monta o dashboard com 6 gráficos e devolve a figura do matplotlib."""
    geral = consultar(conn, "visao_geral")
    fig, eixos = plt.subplots(2, 3, figsize=(18, 10))
    (a1, a2, a3), (a4, a5, a6) = eixos

    # 1. Categorias por gene (barras empilhadas, em %)
    cat = _tabela_categorias(consultar(conn, "categorias_por_gene"), "gene")
    pct = cat.div(cat.sum(axis=1), axis=0) * 100
    pct.plot.barh(stacked=True, ax=a1, color=[CORES_CATEGORIA[c] for c in pct.columns], width=0.6)
    a1.set_title("Classificação clínica por gene (%)")
    a1.set_xlabel("% das variantes")
    a1.set_ylabel("")
    a1.legend(fontsize=8, loc="lower right")

    # 2. Consequências mais comuns entre as patogênicas
    cons = consultar(conn, "consequencias_patogenicas")
    tabela = cons.pivot_table(index="consequencia", columns="gene", values="variantes",
                              aggfunc="sum", fill_value=0)
    tabela = tabela.loc[tabela.sum(axis=1).sort_values().index]
    tabela.plot.barh(ax=a2, color=[CORES_GENE.get(g, "#777") for g in tabela.columns])
    a2.set_title("Consequências mais frequentes\nentre variantes patogênicas")
    a2.set_xlabel("Número de variantes")
    a2.set_ylabel("")
    a2.legend(title="gene", fontsize=8, loc="lower right")

    # 3. Taxa de patogenicidade por consequência
    taxa = consultar(conn, "taxa_patogenicidade_por_consequencia").sort_values("pct_patogenicas")
    a3.barh(taxa["consequencia"], taxa["pct_patogenicas"], color="#C8553D")
    a3.set_xlim(0, 100)
    a3.set_title("Chance de ser patogênica, por consequência\n(entre variantes já classificadas)")
    a3.set_xlabel("% patogênicas")

    # 4. Qualidade da evidência (estrelas)
    est = consultar(conn, "qualidade_evidencia")
    est_tab = est.pivot_table(index="estrelas", columns="categoria", values="variantes",
                              aggfunc="sum", fill_value=0)
    est_tab = est_tab[[c for c in ORDEM_CATEGORIAS if c in est_tab.columns]]
    est_tab.plot.bar(stacked=True, ax=a4, color=[CORES_CATEGORIA[c] for c in est_tab.columns])
    a4.set_title("Nível de revisão (estrelas do ClinVar)")
    a4.set_xlabel("Estrelas (0 = sem critério, 4 = diretriz clínica)")
    a4.set_ylabel("Número de variantes")
    a4.tick_params(axis="x", rotation=0)
    a4.legend(fontsize=8)

    # 5. Evolução anual das avaliações
    evo = consultar(conn, "evolucao_anual")
    evo_tab = _tabela_categorias(evo, "ano")
    for categoria in evo_tab.columns:
        a5.plot(evo_tab.index, evo_tab[categoria], marker="o", ms=3,
                color=CORES_CATEGORIA[categoria], label=categoria)
    a5.set_title("Ano da última avaliação das variantes")
    a5.set_xlabel("Ano")
    a5.set_ylabel("Número de variantes")
    a5.legend(fontsize=8)

    # 6. Hotspots ao longo dos genes
    hot = consultar(conn, "hotspots")
    for gene, dados in hot.groupby("gene"):
        # Descarta blocos isolados nas pontas (0,5% das variantes de cada lado), que
        # correspondem a registros em regiões vizinhas e distorceriam o eixo.
        dados = dados.sort_values("inicio_bloco")
        acumulado = dados["patogenicas"].cumsum() / dados["patogenicas"].sum()
        inicio_faixa = acumulado - dados["patogenicas"] / dados["patogenicas"].sum()
        dados = dados[(acumulado > 0.005) & (inicio_faixa < 0.995)]
        posicao_relativa = (dados["inicio_bloco"] - dados["inicio_bloco"].min()) / 1000
        a6.bar(posicao_relativa, dados["patogenicas"], width=1.8, alpha=0.7,
               color=CORES_GENE.get(gene, "#777"), label=gene)
    a6.set_title("Distribuição das variantes patogênicas\nao longo do gene (blocos de 2 kb)")
    a6.set_xlabel("Distância do início da região (kb)")
    a6.set_ylabel("Variantes patogênicas")
    a6.legend(fontsize=8)

    resumo = "   |   ".join(
        f"{r.gene}: {_br(r.total_variantes)} variantes, {_br(r.pct_patogenicas)}% patogênicas, "
        f"{_br(r.pct_incertas)}% de significado incerto"
        for r in geral.itertuples()
    )
    fig.suptitle("Variantes de BRCA1 e BRCA2 no ClinVar\n" + resumo, fontsize=13)
    for eixo in eixos.flat:
        eixo.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    if salvar_em:
        fig.savefig(salvar_em, dpi=120, bbox_inches="tight")
    return fig
