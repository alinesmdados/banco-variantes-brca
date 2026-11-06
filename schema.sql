-- =============================================================================
-- Banco de variantes genéticas de BRCA1 e BRCA2 (fonte: NCBI ClinVar)
-- Compatível com SQLite 3.
-- =============================================================================

PRAGMA foreign_keys = ON;

DROP VIEW  IF EXISTS vw_variantes;
DROP TABLE IF EXISTS variante_condicao;
DROP TABLE IF EXISTS condicoes;
DROP TABLE IF EXISTS variantes;
DROP TABLE IF EXISTS classificacoes;
DROP TABLE IF EXISTS genes;

-- -----------------------------------------------------------------------------
-- Genes analisados
-- -----------------------------------------------------------------------------
CREATE TABLE genes (
    gene_id     INTEGER PRIMARY KEY,          -- NCBI Gene ID (BRCA1 = 672, BRCA2 = 675)
    simbolo     TEXT    NOT NULL UNIQUE,
    cromossomo  TEXT    NOT NULL,
    descricao   TEXT
);

-- -----------------------------------------------------------------------------
-- Classificação clínica (germinativa) atribuída pelo ClinVar
-- A coluna "categoria" agrupa as descrições originais em 5 grupos de análise.
-- -----------------------------------------------------------------------------
CREATE TABLE classificacoes (
    classificacao_id INTEGER PRIMARY KEY AUTOINCREMENT,
    descricao        TEXT NOT NULL UNIQUE,    -- texto original do ClinVar
    categoria        TEXT NOT NULL CHECK (categoria IN
                     ('Patogênica', 'Benigna', 'Significado incerto', 'Conflitante', 'Outra'))
);

-- -----------------------------------------------------------------------------
-- Variantes
-- -----------------------------------------------------------------------------
CREATE TABLE variantes (
    variante_id        INTEGER PRIMARY KEY,   -- UID do ClinVar
    accession          TEXT    NOT NULL,      -- ex.: VCV000017661
    gene_id            INTEGER NOT NULL REFERENCES genes (gene_id),
    titulo             TEXT,                  -- nomenclatura HGVS
    tipo               TEXT,                  -- ex.: single nucleotide variant, deletion
    consequencia       TEXT,                  -- ex.: missense variant, frameshift variant
    alteracao_proteina TEXT,                  -- ex.: C61G
    posicao_inicio     INTEGER,               -- coordenada GRCh38
    posicao_fim        INTEGER,
    classificacao_id   INTEGER REFERENCES classificacoes (classificacao_id),
    status_revisao     TEXT,
    estrelas           INTEGER CHECK (estrelas BETWEEN 0 AND 4),
    ultima_avaliacao   DATE
);

CREATE INDEX idx_variantes_gene          ON variantes (gene_id);
CREATE INDEX idx_variantes_classificacao ON variantes (classificacao_id);

-- -----------------------------------------------------------------------------
-- Condições (doenças/fenótipos) associadas às variantes: relação N:N
-- -----------------------------------------------------------------------------
CREATE TABLE condicoes (
    condicao_id INTEGER PRIMARY KEY AUTOINCREMENT,
    nome        TEXT NOT NULL UNIQUE
);

CREATE TABLE variante_condicao (
    variante_id INTEGER NOT NULL REFERENCES variantes (variante_id) ON DELETE CASCADE,
    condicao_id INTEGER NOT NULL REFERENCES condicoes (condicao_id),
    PRIMARY KEY (variante_id, condicao_id)
);

-- -----------------------------------------------------------------------------
-- View analítica: uma linha por variante, já com gene e categoria
-- -----------------------------------------------------------------------------
CREATE VIEW vw_variantes AS
SELECT
    v.variante_id,
    v.accession,
    g.simbolo                         AS gene,
    v.titulo,
    v.tipo,
    v.consequencia,
    v.alteracao_proteina,
    v.posicao_inicio,
    c.descricao                       AS classificacao,
    c.categoria,
    v.status_revisao,
    v.estrelas,
    v.ultima_avaliacao,
    CAST(strftime('%Y', v.ultima_avaliacao) AS INTEGER) AS ano_avaliacao
FROM variantes v
JOIN genes g               ON g.gene_id = v.gene_id
LEFT JOIN classificacoes c ON c.classificacao_id = v.classificacao_id;
