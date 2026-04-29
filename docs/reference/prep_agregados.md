# Prepara dados de peneiramento e calcula DMP e DMG

Prepara dados de peneiramento e calcula DMP e DMG

## Usage

``` r
prep_agregados(
  df,
  col_amostra = "amostra_id",
  col_diametro = "diametro_mm",
  col_massa = "massa_g"
)
```

## Arguments

- df:

  Data frame com os dados brutos.

- col_amostra:

  String. Nome da coluna de identificacao da amostra.

- col_diametro:

  String. Nome da coluna de diametro (mm).

- col_massa:

  String. Nome da coluna com a massa (g).

## Value

Uma lista contendo os dados processados e os indices tradicionais.
