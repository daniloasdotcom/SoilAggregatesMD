# Executa analise completa e detalhada de agregacao

Executa analise completa e detalhada de agregacao

## Usage

``` r
analise_completa_dma(
  dados,
  col_amostra = "amostra_id",
  col_diametro = "diametro_mm",
  col_massa = "massa_g"
)
```

## Arguments

- dados:

  Data frame bruto de peneiramento.

- col_amostra:

  String. Nome da coluna de identificacao.

- col_diametro:

  String. Nome da coluna de diametro (mm).

- col_massa:

  String. Nome da coluna com a massa (g).

## Value

Uma lista contendo o resumo consolidado e os detalhes dos ajustes das
equacoes.
