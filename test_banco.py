"""Testes automatizados do ETL e das consultas (execute com: pytest)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import etl_clinvar as etl  # noqa: E402

# Documento no formato atual do esummary do ClinVar (campos relevantes)
DOC_ATUAL = {
    "uid": "17661",
    "accession": "VCV000017661",
    "title": "NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)",
    "obj_type": "single nucleotide variant",
    "protein_change": "C61G, C14G",
    "molecular_consequence_list": ["missense variant"],
    "genes": [{"symbol": "BRCA1", "geneid": "672"}],
    "germline_classification": {
        "description": "Pathogenic",
        "review_status": "reviewed by expert panel",
        "last_evaluated": "2016/10/05 00:00",
        "trait_set": [
            {"trait_name": "Hereditary breast ovarian cancer syndrome"},
            {"trait_name": "not provided"},
        ],
    },
    "variation_set": [{
        "variation_loc": [
            {"assembly_name": "GRCh37", "status": "previous", "start": "41258474", "stop": "41258474"},
            {"assembly_name": "GRCh38", "status": "current", "start": "43106457", "stop": "43106457"},
        ],
    }],
}

# Documento no formato antigo (clinical_significance)
DOC_ANTIGO = {
    "uid": "51000",
    "accession": "VCV000051000",
    "title": "NM_000059.4(BRCA2):c.68-7T>A",
    "genes": [{"symbol": "BRCA2"}],
    "clinical_significance": {
        "description": "Benign",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "last_evaluated": "2021/03/01 00:00",
    },
    "trait_set": [{"trait_name": "Breast-ovarian cancer, familial 2"}],
    "variation_set": [{"variant_type": "single nucleotide variant", "variation_loc": []}],
}


@pytest.mark.parametrize("descricao, esperado", [
    ("Pathogenic", "Patogênica"),
    ("Likely pathogenic", "Patogênica"),
    ("Pathogenic/Likely pathogenic", "Patogênica"),
    ("Benign/Likely benign", "Benigna"),
    ("Uncertain significance", "Significado incerto"),
    ("Conflicting classifications of pathogenicity", "Conflitante"),
    ("not provided", "Outra"),
    (None, "Outra"),
])
def test_categorizar(descricao, esperado):
    assert etl.categorizar(descricao) == esperado


def test_estrelas():
    assert etl.contar_estrelas("practice guideline") == 4
    assert etl.contar_estrelas("reviewed by expert panel") == 3
    assert etl.contar_estrelas("criteria provided, single submitter") == 1
    assert etl.contar_estrelas("no assertion criteria provided") == 0
    assert etl.contar_estrelas(None) == 0


def test_converter_data():
    assert etl.converter_data("2016/10/05 00:00") == "2016-10-05"
    assert etl.converter_data("1/1/2020") is None
    assert etl.converter_data(None) is None


def test_extrair_formato_atual():
    r = etl.extrair_variante(DOC_ATUAL, "BRCA1")
    assert r["variante_id"] == 17661
    assert r["classificacao"] == "Pathogenic"
    assert r["estrelas"] == 3
    assert r["consequencia"] == "missense variant"
    assert r["alteracao_proteina"] == "C61G"
    assert (r["posicao_inicio"], r["posicao_fim"]) == (43106457, 43106457)  # GRCh38
    assert r["ultima_avaliacao"] == "2016-10-05"
    assert r["condicoes"] == ["Hereditary breast ovarian cancer syndrome"]  # "not provided" removido


def test_extrair_formato_antigo():
    r = etl.extrair_variante(DOC_ANTIGO, "BRCA2")
    assert r["classificacao"] == "Benign"
    assert r["tipo"] == "single nucleotide variant"
    assert r["estrelas"] == 2
    assert r["posicao_inicio"] is None
    assert r["condicoes"] == ["Breast-ovarian cancer, familial 2"]


def test_variante_de_outro_gene_e_ignorada():
    assert etl.extrair_variante(DOC_ATUAL, "BRCA2") is None


@pytest.fixture
def banco(tmp_path):
    conn = etl.criar_banco(str(tmp_path / "teste.db"))
    etl.carregar_variantes(conn, "BRCA1", [etl.extrair_variante(DOC_ATUAL, "BRCA1")])
    etl.carregar_variantes(conn, "BRCA2", [etl.extrair_variante(DOC_ANTIGO, "BRCA2")])
    yield conn
    conn.close()


def test_carga_e_integridade(banco):
    assert banco.execute("SELECT COUNT(*) FROM variantes").fetchone()[0] == 2
    assert banco.execute("SELECT COUNT(*) FROM genes").fetchone()[0] == 2
    assert banco.execute("SELECT COUNT(*) FROM condicoes").fetchone()[0] == 2
    assert banco.execute("PRAGMA foreign_key_check").fetchall() == []


def test_carga_e_idempotente(banco):
    n = etl.carregar_variantes(banco, "BRCA1", [etl.extrair_variante(DOC_ATUAL, "BRCA1")])
    assert n == 0
    assert banco.execute("SELECT COUNT(*) FROM variantes").fetchone()[0] == 2


def test_view_e_categoria(banco):
    linhas = dict(banco.execute("SELECT gene, categoria FROM vw_variantes").fetchall())
    assert linhas == {"BRCA1": "Patogênica", "BRCA2": "Benigna"}


def test_todas_as_consultas_executam(banco):
    from dashboard import CONSULTAS
    assert len(CONSULTAS) >= 8
    for nome, sql in CONSULTAS.items():
        banco.execute(sql).fetchall()  # não deve gerar erro
