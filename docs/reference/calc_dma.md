# Calcula o Diametro Medio Ajustado (DMA)

Calcula o Diametro Medio Ajustado (DMA)

## Usage

``` r
calc_dma(
  df_processado,
  col_amostra = "amostra_id",
  col_diametro = "diametro_mm"
)
```

## Arguments

- df_processado:

  Data frame retornado por prep_agregados.

- col_amostra:

  String. Nome da coluna de identificacao.

- col_diametro:

  String. Nome da coluna de diametro (mm).

## Value

Data frame de resumo com os resultados do DMA por amostra.
