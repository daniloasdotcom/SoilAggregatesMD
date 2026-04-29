# Plota a distribuicao de tamanho de agregados

Cria um grafico de barras ilustrando a proporcao de massa retida em cada
classe de diametro para as amostras selecionadas.

## Usage

``` r
plot_distribuicao_agregados(
  df_processado,
  col_amostra = "amostra_id",
  col_diametro = "diametro_mm",
  amostras_selecionadas = NULL
)
```

## Arguments

- df_processado:

  Data frame. O objeto retornado pela funcao
  [`prep_agregados()`](prep_agregados.md).

- col_amostra:

  String. Nome da coluna de identificacao da amostra.

- col_diametro:

  String. Nome da coluna de diametro limite das classes (mm).

- amostras_selecionadas:

  Vetor opcional de caracteres com IDs das amostras a serem plotadas.

## Value

Um objeto ggplot com o painel de distribuicao por amostra.
