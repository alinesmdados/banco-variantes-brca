"""
etl_clinvar.py
==============

Pipeline ETL (Extração, Transformação e Carga) que monta um banco SQLite
com as variantes genéticas de BRCA1 e BRCA2 registradas no NCBI ClinVar.

    Extração      -> API pública E-utilities do NCBI (esearch + esummary, JSON)
    Transformação -> normalização da classificação clínica, estrelas de revisão,
                     datas, coordenadas genômicas (GRCh38) e condições
    Carga         -> banco SQLite modelado em schema.sql

Uso pela linha de comando:
    python etl_clinvar.py --db clinvar_brca.db --genes BRCA1 BRCA2

Autora: Aline Machado
Licença: MIT
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import time
import urllib.parse
import urllib.request
from typing import Iterable, Iterator

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FERRAMENTA = "banco-variantes-brca"
PAUSA_ENTRE_REQUISICOES = 0.4   # o NCBI permite até 3 requisições/s sem chave de API
TAMANHO_LOTE = 400

GENES = {
    "BRCA1": {"gene_id": 672, "cromossomo": "17", "descricao": "BRCA1 DNA repair associated"},
    "BRCA2": {"gene_id": 675, "cromossomo": "13", "descricao": "BRCA2 DNA repair associated"},
}

# Status de revisão do ClinVar -> estrelas (0 a 4), conforme a documentação do ClinVar
ESTRELAS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
}


# ---------------------------------------------------------------------------
# Extração
# ---------------------------------------------------------------------------

def _requisitar(utilidade: str, parametros: dict) -> dict:
    """Faz uma requisição GET às E-utilities e devolve o JSON decodificado."""
    parametros = {**parametros, "retmode": "json", "tool": FERRAMENTA}
    url = f"{EUTILS}/{utilidade}.fcgi?{urllib.parse.urlencode(parametros)}"
    for tentativa in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as resposta:
                dados = json.loads(resposta.read().decode("utf-8"))
            time.sleep(PAUSA_ENTRE_REQUISICOES)
            return dados
        except Exception:  # erros de rede: espera e tenta de novo
            if tentativa == 3:
                raise
            time.sleep(2 * (tentativa + 1))
    raise RuntimeError("não alcançado")


def buscar_variantes(gene: str) -> tuple[int, str, str]:
    """Busca no ClinVar as variantes que afetam apenas o gene informado.

    Retorna (quantidade, WebEnv, query_key): o resultado fica guardado
    no servidor do NCBI ("history server") para ser baixado em lotes.
    """
    dados = _requisitar("esearch", {
        "db": "clinvar",
        "term": f"{gene}[gene] AND single_gene[prop]",
        "usehistory": "y",
        "retmax": 0,
    })
    r = dados["esearchresult"]
    return int(r["count"]), r["webenv"], r["querykey"]


def baixar_resumos(webenv: str, query_key: str, total: int,
                   lote: int = TAMANHO_LOTE) -> Iterator[dict]:
    """Baixa os resumos (esummary) das variantes, em lotes."""
    for inicio in range(0, total, lote):
        dados = _requisitar("esummary", {
            "db": "clinvar", "WebEnv": webenv, "query_key": query_key,
            "retstart": inicio, "retmax": lote,
        })
        resultado = dados.get("result", {})
        for uid in resultado.get("uids", []):
            yield resultado[uid]


# ---------------------------------------------------------------------------
# Transformação
# ---------------------------------------------------------------------------

def categorizar(descricao: str | None) -> str:
    """Agrupa a classificação original do ClinVar em 5 categorias.

    >>> categorizar("Pathogenic/Likely pathogenic")
    'Patogênica'
    >>> categorizar("Likely benign")
    'Benigna'
    >>> categorizar("Conflicting classifications of pathogenicity")
    'Conflitante'
    """
    d = (descricao or "").lower()
    if "conflicting" in d:
        return "Conflitante"
    if "pathogenic" in d and "benign" not in d:
        return "Patogênica"
    if "benign" in d and "pathogenic" not in d:
        return "Benigna"
    if "uncertain" in d:
        return "Significado incerto"
    return "Outra"


def contar_estrelas(status: str | None) -> int:
    """Converte o status de revisão em estrelas (0 a 4).

    >>> contar_estrelas("reviewed by expert panel")
    3
    >>> contar_estrelas("no assertion criteria provided")
    0
    """
    return ESTRELAS.get((status or "").strip().lower(), 0)


def converter_data(texto: str | None) -> str | None:
    """Converte datas do ClinVar ('2023/10/05 00:00') para ISO ('2023-10-05').

    >>> converter_data("2023/10/05 00:00")
    '2023-10-05'
    >>> converter_data("") is None
    True
    """
    m = re.match(r"(\d{4})[/-](\d{2})[/-](\d{2})", texto or "")
    return f"{m.group(1)}-{m.group(2)}-{m.group(3)}" if m else None


def _para_int(valor) -> int | None:
    try:
        return int(str(valor).strip())
    except (TypeError, ValueError):
        return None


def _primeiro(valor) -> str | None:
    """Primeiro item de uma lista ou de um texto separado por vírgulas."""
    if isinstance(valor, list):
        valor = valor[0] if valor else None
    if isinstance(valor, str):
        valor = valor.split(",")[0].strip()
    return valor or None


def extrair_variante(doc: dict, gene: str) -> dict | None:
    """Transforma um documento esummary do ClinVar em um registro do banco.

    Retorna None se a variante não pertencer ao gene informado.
    Aceita tanto o formato atual (germline_classification) quanto o
    antigo (clinical_significance) da API.
    """
    simbolos = {g.get("symbol") for g in doc.get("genes", []) if isinstance(g, dict)}
    if simbolos and gene not in simbolos:
        return None

    bloco = doc.get("germline_classification") or doc.get("clinical_significance") or {}
    variacoes = doc.get("variation_set") or [{}]
    variacao = variacoes[0] if isinstance(variacoes[0], dict) else {}

    inicio = fim = None
    for loc in variacao.get("variation_loc", []) or []:
        if loc.get("assembly_name") == "GRCh38":
            inicio, fim = _para_int(loc.get("start")), _para_int(loc.get("stop"))
            if loc.get("status") == "current":
                break

    condicoes = []
    for trait in bloco.get("trait_set") or doc.get("trait_set") or []:
        nome = (trait.get("trait_name") or "").strip()
        if nome and nome.lower() not in ("not provided", "not specified") and nome not in condicoes:
            condicoes.append(nome)

    return {
        "variante_id": _para_int(doc.get("uid")),
        "accession": doc.get("accession") or f"UID{doc.get('uid')}",
        "gene": gene,
        "titulo": doc.get("title"),
        "tipo": doc.get("obj_type") or variacao.get("variant_type"),
        "consequencia": _primeiro(doc.get("molecular_consequence_list")),
        "alteracao_proteina": _primeiro(doc.get("protein_change")),
        "posicao_inicio": inicio,
        "posicao_fim": fim,
        "classificacao": bloco.get("description") or "not provided",
        "status_revisao": bloco.get("review_status"),
        "estrelas": contar_estrelas(bloco.get("review_status")),
        "ultima_avaliacao": converter_data(bloco.get("last_evaluated")),
        "condicoes": condicoes,
    }


# ---------------------------------------------------------------------------
# Carga
# ---------------------------------------------------------------------------

def criar_banco(caminho_db: str, caminho_schema: str | None = None) -> sqlite3.Connection:
    """Cria (ou recria) o banco a partir do schema.sql."""
    if caminho_schema is None:
        caminho_schema = os.path.join(os.path.dirname(os.path.abspath(__file__)), "schema.sql")
    conn = sqlite3.connect(caminho_db)
    with open(caminho_schema, encoding="utf-8") as f:
        conn.executescript(f.read())
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def _id_classificacao(conn: sqlite3.Connection, descricao: str) -> int:
    conn.execute(
        "INSERT OR IGNORE INTO classificacoes (descricao, categoria) VALUES (?, ?)",
        (descricao, categorizar(descricao)),
    )
    return conn.execute(
        "SELECT classificacao_id FROM classificacoes WHERE descricao = ?", (descricao,)
    ).fetchone()[0]


def _id_condicao(conn: sqlite3.Connection, nome: str) -> int:
    conn.execute("INSERT OR IGNORE INTO condicoes (nome) VALUES (?)", (nome,))
    return conn.execute("SELECT condicao_id FROM condicoes WHERE nome = ?", (nome,)).fetchone()[0]


def carregar_variantes(conn: sqlite3.Connection, gene: str, registros: Iterable[dict]) -> int:
    """Insere o gene e suas variantes no banco. Retorna quantas foram inseridas."""
    info = GENES[gene]
    conn.execute(
        "INSERT OR IGNORE INTO genes (gene_id, simbolo, cromossomo, descricao) VALUES (?, ?, ?, ?)",
        (info["gene_id"], gene, info["cromossomo"], info["descricao"]),
    )
    inseridas = 0
    for r in registros:
        if r is None or r["variante_id"] is None:
            continue
        cursor = conn.execute(
            """INSERT OR IGNORE INTO variantes
               (variante_id, accession, gene_id, titulo, tipo, consequencia, alteracao_proteina,
                posicao_inicio, posicao_fim, classificacao_id, status_revisao, estrelas,
                ultima_avaliacao)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (r["variante_id"], r["accession"], info["gene_id"], r["titulo"], r["tipo"],
             r["consequencia"], r["alteracao_proteina"], r["posicao_inicio"], r["posicao_fim"],
             _id_classificacao(conn, r["classificacao"]), r["status_revisao"], r["estrelas"],
             r["ultima_avaliacao"]),
        )
        if cursor.rowcount:
            inseridas += 1
            for nome in r["condicoes"]:
                conn.execute(
                    "INSERT OR IGNORE INTO variante_condicao (variante_id, condicao_id) VALUES (?, ?)",
                    (r["variante_id"], _id_condicao(conn, nome)),
                )
    conn.commit()
    return inseridas


# ---------------------------------------------------------------------------
# Pipeline completo
# ---------------------------------------------------------------------------

def construir_banco(caminho_db: str = "clinvar_brca.db",
                    genes: Iterable[str] = ("BRCA1", "BRCA2"),
                    verbose: bool = True) -> sqlite3.Connection:
    """Executa o ETL completo e devolve a conexão com o banco pronto."""
    conn = criar_banco(caminho_db)
    for gene in genes:
        total, webenv, query_key = buscar_variantes(gene)
        if verbose:
            print(f"{gene}: {total:,} variantes encontradas no ClinVar. Baixando...".replace(",", "."))
        registros = (extrair_variante(doc, gene) for doc in baixar_resumos(webenv, query_key, total))
        n = carregar_variantes(conn, gene, registros)
        if verbose:
            print(f"{gene}: {n:,} variantes carregadas no banco.".replace(",", "."))
    return conn


def main() -> None:
    parser = argparse.ArgumentParser(description="Monta o banco de variantes BRCA a partir do ClinVar.")
    parser.add_argument("--db", default="clinvar_brca.db", help="arquivo SQLite de saída")
    parser.add_argument("--genes", nargs="+", default=["BRCA1", "BRCA2"], choices=sorted(GENES))
    args = parser.parse_args()
    construir_banco(args.db, args.genes).close()


if __name__ == "__main__":
    main()
