# SoilAggregatesMD 0.2.0

* **Expansão de Modelos para DMA:** O cálculo do Diâmetro Médio Ajustado (`calc_dma()`) agora avalia automaticamente 7 modelos matemáticos. Foram adicionados 4 novos modelos paramétricos globais contínuos (Weibull, Logística, Log-Normal e Gompertz) para representar melhor amostras com assimetrias ou distribuições em "S", somando-se às 3 equações empíricas originais.
* **Diagnóstico Visual Aprimorado:** A função `plot_diagnostico_amostra()` foi atualizada para renderizar uma grade comparativa completa, exibindo os pontos observados, as curvas teóricas, os parâmetros e o R² para todos os 7 modelos simultaneamente.
* **Relatório Detalhado:** A função `analise_completa_dma()` agora consolida e seleciona o modelo de melhor ajuste dentro do novo leque estendido de equações, extraindo o DMA mais preciso para cada amostra.
* **Rigor Metodológico Reforçado:** O pacote consolida a adesão aos pressupostos teóricos de van Lier & Albuquerque (1997), rejeitando técnicas de sobreajuste local (como *Splines*) e mantendo o foco em modelos que traduzem uma distribuição física global da amostra. O alerta para quando o R² máximo for inferior a 0,99 foi mantido.

# SoilAggregatesMD 0.1.0

* Primeira versão oficial do pacote.
* Adicionada a função `prep_agregados()` para processamento inicial de dados de peneiramento.
* Implementado o cálculo do DMA (Diâmetro Médio Ajustado) através de três modelos de regressão (Equações 3, 4 e 5 de van Lier & Albuquerque, 1997).
* Adicionada a função `analise_completa_dma()` para relatórios automatizados.
* Incluídas ferramentas de diagnóstico visual com `plot_diagnostico_amostra()`.
* Implementada a recomendação de segurança de R² >= 0,99.
