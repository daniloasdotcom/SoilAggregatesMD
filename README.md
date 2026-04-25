# SoilAggregatesMD

O objetivo do **SoilAggregatesMD** é automatizar o cálculo e a análise da distribuição de agregados do solo, com foco especial no **Diâmetro Médio Ajustado (DMA)**, conforme proposto por van Lier & Albuquerque (1997).

Este pacote permite processar grandes volumes de dados de peneiramento, testar modelos de regressão não-linear e gerar diagnósticos visuais de precisão para pesquisadores da Ciência do Solo.

## Instalação

Você pode instalar a versão de desenvolvimento do **SoilAggregatesMD** diretamente do GitHub com o seguinte comando no R:

```r
# install.packages("remotes")
remotes::install_github("daniloasdotcom/SoilAggregatesMD")
```

## Funcionalidades Principais

* **Cálculo de Índices Tradicionais:** DMP (Diâmetro Médio Ponderado) e DMG (Diâmetro Médio Geométrico).
* **Modelagem de DMA:** Ajuste automático das Equações 3, 4 e 5 do artigo original.
* **Validação Científica:** Mensagens de aviso automáticas para ajustes com $R^2 < 0,99$.
* **Diagnóstico Visual:** Gráficos comparativos entre os três modelos para cada amostra individual.
* **Exportação Facilitada:** Estrutura compatível com exportação para Excel (.xlsx).

## Exemplo de Uso

Este é um exemplo básico de como processar seus dados usando a função de análise completa:

```r
library(SoilAggregatesMD)
library(readxl)

# 1. Carregar seus dados (formato longo)
dados <- read_excel("seus_dados_solo.xlsx")

# 2. Executar a análise completa (Índices + DMA + Detalhes dos Modelos)
analise <- analise_completa_dma(dados)

# 3. Visualizar o resumo dos índices calculados
print(analise$resumo)

# 4. Gerar gráfico de diagnóstico para uma amostra específica
plot_diagnostico_amostra(dados, amostra_id_alvo = "Amostra_A")
```

## Recomendação Importante

Segundo os autores originais (van Lier & Albuquerque, 1997), o método do Diâmetro Médio Ajustado (DMA) deve ser aplicado sem ressalvas apenas quando pelo menos uma das equações de regressão apresentar um coeficiente de determinação ($R^2$) maior ou igual a 0,99. O pacote emite alertas automáticos caso esse critério não seja atingido.

## Referência Bibliográfica

Van Lier, Q. D. J., & Albuquerque, J. A. (1997). Novo método para calcular o diâmetro médio de agregados de solos. Revista Brasileira De Ciência Do Solo, 21(4), 699–705. https://doi.org/10.1590/S0100-06831997000400022

---
