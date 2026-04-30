# SoilAggregatesMD 0.5.0

* **Arquitetura Orientada a Objetos (Classes S3):** O pacote foi reestruturado para utilizar o sistema de classes S3 do R. A função principal `analise_completa_dma()` agora gera e retorna um objeto da classe `agregados_dma`.
* **Novos Métodos Genéricos nativos:** Adicionadas as funções S3 `summary.agregados_dma()` e `plot.agregados_dma()`. Agora o usuário pode inspecionar a tabela formatada e gerar os gráficos das curvas diretamente com os comandos clássicos `summary(resultado)` e `plot(resultado)`.
* **Robustez Matemática Avançada (Levenberg-Marquardt):** Substituição do otimizador padrão `stats::nls` pela função `minpack.lm::nlsLM` em todas as rotinas de modelagem. O algoritmo reduz significativamente falhas de convergência e singularidades matriciais.
* **Critério de Seleção Híbrido (RMSE + R²):** A mecânica de seleção do "Modelo Vencedor" foi aprimorada e agora calcula a Raiz do Erro Quadrático Médio (RMSE). O pacote elege o modelo que possuir o menor RMSE (maior precisão física no traçado), desde que este respeite o filtro rigoroso de $R^2 \ge 0,99$ exigido pela metodologia de van Lier & Albuquerque (1997)[cite: 1].
* **Novo Índice (Fração > 2 mm):** A função `prep_agregados()` passou a calcular e retornar automaticamente a porcentagem retida maior que 2 mm (`Fracao_Maior_2mm_pct`), um índice simples de agregação que apresenta alta correlação com o DMP[cite: 1].
* **Validação de Dados (Sanity Checks):** Inseridas defesas e validações iniciais em `prep_agregados()` para verificar a estrutura do *data.frame*, a existência das colunas informadas e bloquear valores de massa negativos ou em branco (`NA`), retornando mensagens de erro claras ao usuário.
* **Nova Dependência:** Inclusão do pacote `minpack.lm` na seção *Imports*.

# SoilAggregatesMD 0.4.0

* **Visualização da Distribuição Original:** Adicionada a função `plot_distribuicao_agregados()`, que gera gráficos de barras facetados ilustrando a proporção percentual de massa retida em cada classe de peneira. O gráfico ordena automaticamente as classes do maior para o menor diâmetro, auxiliando na inspeção inicial de qualidade e detecção de assimetrias antes do cálculo dos modelos.

# SoilAggregatesMD 0.3.0

* **Nova Função de Exportação:** Adicionada a função `exportar_analise_xlsx()`, que permite exportar os resultados da `analise_completa_dma()` diretamente para um arquivo Excel (`.xlsx`). O arquivo gerado organiza automaticamente o resumo geral e os detalhes de ajuste dos modelos em abas separadas.
* **Nova Dependência:** Inclusão do pacote `writexl` na seção *Imports* do projeto, garantindo uma exportação nativa, rápida e sem necessidade de instalação do Java na máquina do usuário.
* **Correção de Bug (Gráficos):** Corrigido um erro de "comprimento zero" na função `plot_dma()`. A função `analise_completa_dma()` foi atualizada para repassar corretamente as colunas `Modelo_Vencedor`, `Param_a` e `Param_b` na tabela de resumo, garantindo a renderização perfeita da curva teórica.

# SoilAggregatesMD 0.2.0

* **Expansão de Modelos para DMA:** O cálculo do Diâmetro Médio Ajustado (`calc_dma()`) agora avalia automaticamente 7 modelos matemáticos. Foram adicionados 4 novos modelos paramétricos globais contínuos (Weibull, Logística, Log-Normal e Gompertz) para representar melhor amostras com assimetrias ou distribuições em "S", somando-se às 3 equações empíricas originais.
* **Diagnóstico Visual Aprimorado:** A função `plot_diagnostico_amostra()` foi atualizada para renderizar uma grade comparativa completa, exibindo os pontos observados, as curvas teóricas, os parâmetros e o R² para todos os 7 modelos simultaneamente.
* **Relatório Detalhado:** A função `analise_completa_dma()` agora consolida e seleciona o modelo de melhor ajuste dentro do novo leque estendido de equações, extraindo o DMA mais preciso para cada amostra.
* **Rigor Metodológico Reforçado:** O pacote consolida a adesão aos pressupostos teóricos de van Lier & Albuquerque (1997), rejeitando técnicas de sobreajuste local (como *Splines*) e mantendo o foco em modelos que traduzem uma distribuição física global da amostra[cite: 1]. O alerta para quando o R² máximo for inferior a 0,99 foi mantido[cite: 1].

# SoilAggregatesMD 0.1.0

* Primeira versão oficial do pacote.
* Adicionada a função `prep_agregados()` para processamento inicial de dados de peneiramento.
* Implementado o cálculo do DMA (Diâmetro Médio Ajustado) através de três modelos de regressão (Equações 3, 4 e 5 de van Lier & Albuquerque, 1997)[cite: 1].
* Adicionada a função `analise_completa_dma()` para relatórios automatizados.
* Incluídas ferramentas de diagnóstico visual com `plot_diagnostico_amostra()`.
* Implementada a recomendação de segurança de R² >= 0,99[cite: 1].
