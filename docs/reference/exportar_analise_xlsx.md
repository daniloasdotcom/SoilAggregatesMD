# Exporta os resultados da analise de agregacao para Excel

Exporta os resultados da analise de agregacao para Excel

## Usage

``` r
exportar_analise_xlsx(lista_analise, caminho_arquivo = "resultados_dma.xlsx")
```

## Arguments

- lista_analise:

  Lista. O objeto retornado por
  [`analise_completa_dma()`](analise_completa_dma.md).

- caminho_arquivo:

  String. O nome ou caminho completo do arquivo a ser criado.

## Value

A funcao nao retorna um objeto no R, mas salva o arquivo no disco.
