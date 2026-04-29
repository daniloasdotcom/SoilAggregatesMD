# Plota as curvas de retencao e o ajuste do DMA

Plota as curvas de retencao e o ajuste do DMA

## Usage

``` r
plot_dma(
  df_processado,
  df_dma,
  col_amostra = "amostra_id",
  col_diametro = "diametro_mm",
  amostras_selecionadas = NULL
)
```

## Arguments

- df_processado:

  Data frame retornado por [`prep_agregados()`](prep_agregados.md).

- df_dma:

  Data frame de resumo retornado por
  [`analise_completa_dma()`](analise_completa_dma.md).

- col_amostra:

  String. Nome da coluna de identificacao.

- col_diametro:

  String. Nome da coluna de diametro (mm).

- amostras_selecionadas:

  Vetor opcional de caracteres. IDs das amostras a serem plotadas.

## Value

Um objeto ggplot2 com os graficos.
