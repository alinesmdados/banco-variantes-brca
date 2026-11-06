-- =============================================================================
-- Consultas analíticas do banco de variantes BRCA1/BRCA2
-- Cada consulta começa com "-- nome:" e é carregada automaticamente pelo
-- dashboard.py (função carregar_consultas).
-- =============================================================================

-- nome: visao_geral
-- Quantas variantes existem por gene e qual a proporção de patogênicas?
SELECT
    gene,
    COUNT(*)                                                       AS total_variantes,
    SUM(categoria = 'Patogênica')                                  AS patogenicas,
    ROUND(100.0 * SUM(categoria = 'Patogênica') / COUNT(*), 1)     AS pct_patogenicas,
    SUM(categoria = 'Significado incerto')                         AS incertas,
    ROUND(100.0 * SUM(categoria = 'Significado incerto') / COUNT(*), 1) AS pct_incertas
FROM vw_variantes
GROUP BY gene
ORDER BY gene;

-- nome: categorias_por_gene
-- Distribuição das categorias de classificação clínica em cada gene.
SELECT gene, categoria, COUNT(*) AS variantes
FROM vw_variantes
GROUP BY gene, categoria
ORDER BY gene, variantes DESC;

-- nome: consequencias_patogenicas
-- Quais tipos de consequência molecular mais aparecem entre as variantes patogênicas?
-- Usa função de janela (RANK) para ordenar dentro de cada gene.
SELECT *
FROM (
    SELECT
        gene,
        COALESCE(consequencia, 'não informada') AS consequencia,
        COUNT(*)                                AS variantes,
        RANK() OVER (PARTITION BY gene ORDER BY COUNT(*) DESC) AS ranking
    FROM vw_variantes
    WHERE categoria = 'Patogênica'
    GROUP BY gene, COALESCE(consequencia, 'não informada')
)
WHERE ranking <= 6
ORDER BY gene, ranking;

-- nome: taxa_patogenicidade_por_consequencia
-- Entre as variantes já classificadas (patogênica ou benigna), qual a chance
-- de cada tipo de consequência ser patogênico?
SELECT
    COALESCE(consequencia, 'não informada')                       AS consequencia,
    COUNT(*)                                                      AS classificadas,
    ROUND(100.0 * SUM(categoria = 'Patogênica') / COUNT(*), 1)    AS pct_patogenicas
FROM vw_variantes
WHERE categoria IN ('Patogênica', 'Benigna')
GROUP BY COALESCE(consequencia, 'não informada')
HAVING COUNT(*) >= 20
ORDER BY pct_patogenicas DESC;

-- nome: qualidade_evidencia
-- Nível de revisão (estrelas do ClinVar) por categoria: quão confiáveis são as classificações?
SELECT
    categoria,
    estrelas,
    COUNT(*) AS variantes
FROM vw_variantes
GROUP BY categoria, estrelas
ORDER BY categoria, estrelas;

-- nome: evolucao_anual
-- Quantas variantes tiveram a última avaliação em cada ano, por categoria.
SELECT ano_avaliacao AS ano, categoria, COUNT(*) AS variantes
FROM vw_variantes
WHERE ano_avaliacao IS NOT NULL
GROUP BY ano_avaliacao, categoria
ORDER BY ano_avaliacao;

-- nome: hotspots
-- Regiões de cada gene (blocos de 2.000 pb) com mais variantes patogênicas.
SELECT
    gene,
    (posicao_inicio / 2000) * 2000          AS inicio_bloco,
    COUNT(*)                                AS patogenicas
FROM vw_variantes
WHERE categoria = 'Patogênica' AND posicao_inicio IS NOT NULL
GROUP BY gene, inicio_bloco
ORDER BY gene, inicio_bloco;

-- nome: condicoes_mais_citadas
-- Condições (doenças) mais associadas às variantes patogênicas.
SELECT
    co.nome                       AS condicao,
    COUNT(DISTINCT vc.variante_id) AS variantes_patogenicas
FROM variante_condicao vc
JOIN condicoes   co ON co.condicao_id = vc.condicao_id
JOIN vw_variantes v ON v.variante_id = vc.variante_id
WHERE v.categoria = 'Patogênica'
GROUP BY co.nome
ORDER BY variantes_patogenicas DESC
LIMIT 10;

-- nome: variantes_expert_panel
-- Exemplos de variantes patogênicas revisadas por painel de especialistas (3+ estrelas).
SELECT gene, accession, alteracao_proteina, consequencia, classificacao, ultima_avaliacao
FROM vw_variantes
WHERE categoria = 'Patogênica' AND estrelas >= 3 AND alteracao_proteina IS NOT NULL
ORDER BY ultima_avaliacao DESC
LIMIT 15;
