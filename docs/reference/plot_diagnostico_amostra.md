# Grafico de diagnostico comparativo para uma amostra

Grafico de diagnostico comparativo para uma amostra

## Usage

``` r
plot_diagnostico_amostra(
  dados,
  amostra_id_alvo,
  col_amostra = "amostra_id",
  col_diametro = "diametro_mm",
  col_massa = "massa_g"
)
```

## Arguments

- dados:

  Data frame bruto.

- amostra_id_alvo:

  String. ID da amostra que sera avaliada.

- col_amostra:

  String. Nome da coluna de identificacao.

- col_diametro:

  String. Nome da coluna de diametro (mm).

- col_massa:

  String. Nome da coluna com a massa (g).

## Value

Objeto ggplot com os 7 modelos testados para a mesma amostra.
