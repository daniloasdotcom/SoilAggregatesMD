# SoilAggregatesMD 0.1.0

* Primeira versão oficial do pacote.
* Adicionada a função `prep_agregados()` para processamento inicial de dados de peneiramento.
* Implementado o cálculo do DMA (Diâmetro Médio Ajustado) através de três modelos de regressão (Equações 3, 4 e 5 de van Lier & Albuquerque, 1997).
* Adicionada a função `analise_completa_dma()` para relatórios automatizados.
* Incluídas ferramentas de diagnóstico visual com `plot_diagnostico_amostra()`.
* Implementada a recomendação de segurança de $R^2 \ge 0,99$.
