#' Prepara dados de peneiramento e calcula DMP e DMG
#'
#' Esta função recebe um data frame no formato longo contendo dados de
#' agregados de solo e calcula as frações de massa, o diâmetro aritmético
#' e os índices tradicionais (DMP e DMG) para cada amostra.
#'
#' @param df Data frame contendo os dados de peneiramento.
#' @param col_amostra String. Nome da coluna de identificação da amostra.
#' @param col_diametro String. Nome da coluna de diâmetro das peneiras (mm).
#' @param col_massa String. Nome da coluna de massa retida (g).
#'
#' @return Uma lista contendo dois data frames: \code{dados_processados}
#' (com as frações calculadas) e \code{indices} (com DMP e DMG por amostra).
#' @export
#'
prep_agregados <- function(df, col_amostra = "amostra_id", col_diametro = "diametro_mm", col_massa = "massa_g") {

  # 1. Cálculos de fração e F acumulada agrupados por amostra
  tabela_base <- df |>
    dplyr::group_by(.data[[col_amostra]]) |>
    dplyr::arrange(dplyr::desc(.data[[col_diametro]]), .by_group = TRUE) |>
    dplyr::mutate(
      massa_total = sum(.data[[col_massa]], na.rm = TRUE),
      fracao_massa = .data[[col_massa]] / massa_total,
      diametro_inf = dplyr::lead(.data[[col_diametro]], default = 0),
      D_aritmetico = (.data[[col_diametro]] + diametro_inf) / 2
    ) |>
    dplyr::arrange(.data[[col_diametro]], .by_group = TRUE) |>
    dplyr::mutate(F_acumulada = cumsum(fracao_massa)) |>
    dplyr::arrange(dplyr::desc(.data[[col_diametro]]), .by_group = TRUE) |>
    dplyr::mutate(
      F_acumulada_inf = dplyr::lead(F_acumulada, default = 0),
      F_media = (F_acumulada + F_acumulada_inf) / 2
    ) |>
    dplyr::ungroup()

  # 2. Cálculo dos índices tradicionais (DMP e DMG) por amostra
  indices_tradicionais <- tabela_base |>
    dplyr::group_by(.data[[col_amostra]]) |>
    dplyr::summarise(
      DMP = sum(fracao_massa * D_aritmetico, na.rm = TRUE),
      DMG = exp(sum(fracao_massa * log(D_aritmetico), na.rm = TRUE)),
      .groups = "drop"
    )

  return(list(
    dados_processados = tabela_base,
    indices = indices_tradicionais
  ))
}


#' Calcula o Diâmetro Médio Ajustado (DMA)
#'
#' Esta função recebe os dados processados e aplica regressões não-lineares
#' para encontrar o melhor ajuste (maior R2) entre as equações propostas por
#' van Lier & Albuquerque (1997). Em seguida, calcula o DMA.
#'
#' @note Segundo van Lier & Albuquerque (1997), o método proposto para o cálculo
#' do DMA deve ser aplicado sem ressalvas apenas quando o coeficiente de
#' determinação (R2) do ajuste da equação for maior ou igual a 0,99.
#'
#' @param df_processado Data frame contendo o elemento `dados_processados`
#' gerado pela função \code{\link{prep_agregados}}.
#' @param col_amostra String. Nome da coluna de identificação da amostra.
#' @param col_diametro String. Nome da coluna de diâmetro das peneiras (mm).
#'
#' @return Um data frame contendo o DMA, o modelo selecionado e os
#' parâmetros (R2, a, b) para cada amostra.
#' @export
#'
calc_dma <- function(df_processado, col_amostra = "amostra_id", col_diametro = "diametro_mm") {

  # Função interna auxiliar para calcular o R2 do modelo nls
  calcular_r2 <- function(observado, predito) {
    sq_tot <- sum((observado - mean(observado, na.rm = TRUE))^2, na.rm = TRUE)
    sq_res <- sum((observado - predito)^2, na.rm = TRUE)
    return(1 - (sq_res / sq_tot))
  }

  # Divide os dados em uma lista, separando por amostra
  amostras_lista <- split(df_processado, df_processado[[col_amostra]])

  # Loop sobre cada amostra usando lapply (robusto e não precisa de pacotes extras)
  resultados_dma <- lapply(amostras_lista, function(dados_amostra) {

    nome_amostra <- unique(dados_amostra[[col_amostra]])

    # Extrai os vetores para evitar problemas de escopo dentro do nls()
    D_vetor <- dados_amostra[[col_diametro]]
    F_vetor <- dados_amostra$F_acumulada

    modelos <- list()

    # Eq 3: F = 1 / (a + b/D)
    tryCatch({
      mod1 <- stats::nls(F_vetor ~ 1 / (a + (b / D_vetor)), start = list(a = 1, b = 1))
      modelos[["Eq3"]] <- list(nome = "Eq3", a = stats::coef(mod1)["a"], b = stats::coef(mod1)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod1)))
    }, error = function(e) NULL)

    # Eq 4: F = 1 / (a + b/D^0.5)
    tryCatch({
      mod2 <- stats::nls(F_vetor ~ 1 / (a + (b / sqrt(D_vetor))), start = list(a = 1, b = 1))
      modelos[["Eq4"]] <- list(nome = "Eq4", a = stats::coef(mod2)["a"], b = stats::coef(mod2)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod2)))
    }, error = function(e) NULL)

    # Eq 5: F = a * D^b
    tryCatch({
      mod3 <- stats::nls(F_vetor ~ a * (D_vetor^b), start = list(a = 1, b = 1))
      modelos[["Eq5"]] <- list(nome = "Eq5", a = stats::coef(mod3)["a"], b = stats::coef(mod3)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod3)))
    }, error = function(e) NULL)

    # Filtra apenas os modelos que convergiram com sucesso
    modelos <- Filter(Negate(is.null), modelos)

    # Tratativa caso nenhum modelo convirja
    if (length(modelos) == 0) {
      return(data.frame(amostra_id = nome_amostra, DMA = NA, Modelo = NA, R2 = NA, a = NA, b = NA))
    }

    # Seleciona o melhor modelo (Maior R2)
    melhor <- modelos[[which.max(sapply(modelos, function(x) x$r2))]]

    # --- CÁLCULO DO DIÂMETRO AJUSTADO ---
    soma_F <- dados_amostra$F_acumulada + dados_amostra$F_acumulada_inf

    if (melhor$nome == "Eq3") {
      D_ajustado <- (melhor$b * soma_F) / (2 - (melhor$a * soma_F))
    } else if (melhor$nome == "Eq4") {
      D_ajustado <- ((melhor$b * soma_F) / (2 - (melhor$a * soma_F)))^2
    } else {
      D_ajustado <- (soma_F / (2 * melhor$a))^(1 / melhor$b)
    }

    # Ponderação pela massa para chegar no DMA
    # Se der infinito ou NaN por divisão por zero, substituímos por 0
    D_ajustado[!is.finite(D_ajustado)] <- 0
    D_ponderado <- dados_amostra$fracao_massa * D_ajustado
    DMA_final <- sum(D_ponderado, na.rm = TRUE)

    # Monta a linha de resultado para esta amostra
    df_resumo <- data.frame(
      amostra_id = nome_amostra,
      DMA = DMA_final,
      Modelo_Vencedor = melhor$nome,
      R2 = melhor$r2,
      Param_a = melhor$a,
      Param_b = melhor$b
    )

    # Renomeia a coluna amostra_id para o nome original que o usuário passou
    colnames(df_resumo)[1] <- col_amostra

    return(df_resumo)
  })

  # Empilha os resultados de todas as amostras em uma única tabela
  tabela_final <- do.call(rbind, resultados_dma)
  rownames(tabela_final) <- NULL

  # Verificação da regra de van Lier & Albuquerque (1997)
  if (any(tabela_final$R2 < 0.99, na.rm = TRUE)) {
    warning("Atenção: Pelo menos uma amostra apresentou R2 máximo inferior a 0,99. ",
            "van Lier & Albuquerque (1997) recomendam o uso do DMA apenas quando R2 >= 0,99.")
  }

  return(tabela_final)
}


#' Plota as curvas de retenção e o ajuste do DMA
#'
#' Gera um gráfico facetado por amostra, mostrando os pontos observados,
#' a curva teórica do modelo vencedor e uma linha indicando o DMA.
#'
#' @note Segundo van Lier & Albuquerque (1997), o método proposto para o cálculo
#' do DMA deve ser aplicado sem ressalvas apenas quando o coeficiente de
#' determinação (R2) do ajuste da equação for maior ou igual a 0,99.
#'
#' @param df_processado Data frame gerado pela função prep_agregados.
#' @param df_dma Data frame gerado pela função calc_dma.
#' @param col_amostra String. Nome da coluna de identificação da amostra.
#' @param col_diametro String. Nome da coluna de diâmetro.
#'
#' @return Um objeto ggplot.
#' @export
#'
plot_dma <- function(df_processado, df_dma, col_amostra = "amostra_id", col_diametro = "diametro_mm") {

  # Função interna para gerar pontos da linha suave da equação
  gerar_curva <- function(amostra, mod, param_a, param_b, d_max) {
    # Cria 100 pontos entre o mínimo e o máximo diâmetro
    d_seq <- seq(0.01, d_max, length.out = 100)

    if (mod == "Eq3") {
      f_seq <- 1 / (param_a + (param_b / d_seq))
    } else if (mod == "Eq4") {
      f_seq <- 1 / (param_a + (param_b / sqrt(d_seq)))
    } else {
      f_seq <- param_a * (d_seq^param_b)
    }

    df <- data.frame(amostra_id = amostra, diametro_mm = d_seq, F_teorica = f_seq)
    # Garante que o nome da coluna da amostra seja dinâmico
    colnames(df)[1] <- col_amostra
    colnames(df)[2] <- col_diametro
    return(df)
  }

  # Aplica a função para cada linha (amostra) da tabela de resultados
  d_maximo <- max(df_processado[[col_diametro]], na.rm = TRUE)

  lista_curvas <- lapply(1:nrow(df_dma), function(i) {
    gerar_curva(
      amostra = df_dma[[col_amostra]][i],
      mod = df_dma$Modelo_Vencedor[i],
      param_a = df_dma$Param_a[i],
      param_b = df_dma$Param_b[i],
      d_max = d_maximo
    )
  })

  # Empilha todas as curvas teóricas geradas
  df_curvas <- do.call(rbind, lista_curvas)

  # Constrói o gráfico com ggplot2
  p <- ggplot2::ggplot() +
    # A linha suave do modelo teórico (azul)
    ggplot2::geom_line(data = df_curvas, ggplot2::aes(x = .data[[col_diametro]], y = F_teorica), color = "blue", linewidth = 0.8) +
    # Os pontos reais observados (preto)
    ggplot2::geom_point(data = df_processado, ggplot2::aes(x = .data[[col_diametro]], y = F_acumulada), size = 2) +
    # Linha vertical indicando exatamente onde o DMA caiu (vermelho)
    ggplot2::geom_vline(data = df_dma, ggplot2::aes(xintercept = DMA), color = "darkred", linetype = "dashed") +
    # Divide o gráfico em vários quadros (um por amostra)
    ggplot2::facet_wrap(stats::as.formula(paste("~", col_amostra))) +
    # Estética final
    ggplot2::labs(
      title = "Ajuste do Diâmetro Médio de Agregados (DMA)",
      subtitle = "Linha azul: Modelo Vencedor | Linha tracejada vermelha: DMA",
      x = "Diâmetro (mm)",
      y = "Fração Acumulada (g/g)",
      caption = "Nota: van Lier & Albuquerque (1997) recomendam o uso do DMA apenas se R² >= 0,99."
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}

#' Executa análise completa e detalhada de agregação
#'
#' Esta função automatiza todo o fluxo: prepara os dados, calcula os índices
#' tradicionais e o DMA, e fornece um relatório detalhado do desempenho de
#' todas as três equações de ajuste para cada amostra, incluindo o valor do
#' DMA que cada uma resultaria.
#'
#' @note Segundo van Lier & Albuquerque (1997), o método proposto para o cálculo
#' do DMA deve ser aplicado sem ressalvas apenas quando o coeficiente de
#' determinação (R2) do ajuste da equação for maior ou igual a 0,99.
#'
#' @param dados Data frame bruto (formato longo) com os dados de peneiramento.
#' @param col_amostra String. Nome da coluna de amostra.
#' @param col_diametro String. Nome da coluna de diâmetro (mm).
#' @param col_massa String. Nome da coluna de massa (g).
#'
#' @return Uma lista contendo:
#' \itemize{
#'   \item \code{resumo}: Tabela com DMP, DMG, DMA, o melhor modelo e o status de recomendação para cada amostra.
#'   \item \code{detalhes_equacoes}: Tabela detalhada com R2, parâmetros a, b e o DMA resultante de TODAS as equações testadas.
#' }
#' @export
#'
analise_completa_dma <- function(dados, col_amostra = "amostra_id", col_diametro = "diametro_mm", col_massa = "massa_g") {

  # 1. Preparação interna
  prep <- prep_agregados(dados, col_amostra, col_diametro, col_massa)
  df_proc <- prep$dados_processados
  indices_trad <- prep$indices

  # Função auxiliar para R2
  calc_r2 <- function(obs, pred) {
    sq_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    sq_res <- sum((obs - pred)^2, na.rm = TRUE)
    return(1 - (sq_res / sq_tot))
  }

  # Função auxiliar para calcular DMA específico de uma equação
  calc_dma_especifico <- function(dados_st, eq_nome, a, b) {
    soma_F <- dados_st$F_acumulada + dados_st$F_acumulada_inf
    if (eq_nome == "Eq3") {
      D_aj <- (b * soma_F) / (2 - (a * soma_F))
    } else if (eq_nome == "Eq4") {
      D_aj <- ((b * soma_F) / (2 - (a * soma_F)))^2
    } else {
      D_aj <- (soma_F / (2 * a))^(1 / b)
    }
    D_aj[!is.finite(D_aj)] <- 0
    return(sum(dados_st$fracao_massa * D_aj, na.rm = TRUE))
  }

  amostras_lista <- split(df_proc, df_proc[[col_amostra]])

  lista_resumo_dma <- list()
  lista_detalhes_modelos <- list()

  for (i in seq_along(amostras_lista)) {
    dados_st <- amostras_lista[[i]]
    id <- unique(dados_st[[col_amostra]])
    D <- dados_st[[col_diametro]]
    F_obs <- dados_st$F_acumulada

    resultados_st <- data.frame()

    # --- Testar Eq 3 ---
    tryCatch({
      m <- stats::nls(F_obs ~ 1 / (a + (b / D)), start = list(a = 1, b = 1))
      a_val <- stats::coef(m)["a"]; b_val <- stats::coef(m)["b"]
      dma_eq <- calc_dma_especifico(dados_st, "Eq3", a_val, b_val)
      resultados_st <- rbind(resultados_st, data.frame(
        amostra_id = id, Equacao = "Eq3", R2 = calc_r2(F_obs, stats::predict(m)),
        DMA_calculado = dma_eq, Param_a = a_val, Param_b = b_val
      ))
    }, error = function(e) {
      resultados_st <<- rbind(resultados_st, data.frame(amostra_id = id, Equacao = "Eq3", R2 = NA, DMA_calculado = NA, Param_a = NA, Param_b = NA))
    })

    # --- Testar Eq 4 ---
    tryCatch({
      m <- stats::nls(F_obs ~ 1 / (a + (b / sqrt(D))), start = list(a = 1, b = 1))
      a_val <- stats::coef(m)["a"]; b_val <- stats::coef(m)["b"]
      dma_eq <- calc_dma_especifico(dados_st, "Eq4", a_val, b_val)
      resultados_st <- rbind(resultados_st, data.frame(
        amostra_id = id, Equacao = "Eq4", R2 = calc_r2(F_obs, stats::predict(m)),
        DMA_calculado = dma_eq, Param_a = a_val, Param_b = b_val
      ))
    }, error = function(e) {
      resultados_st <<- rbind(resultados_st, data.frame(amostra_id = id, Equacao = "Eq4", R2 = NA, DMA_calculado = NA, Param_a = NA, Param_b = NA))
    })

    # --- Testar Eq 5 ---
    tryCatch({
      m <- stats::nls(F_obs ~ a * (D^b), start = list(a = 1, b = 1))
      a_val <- stats::coef(m)["a"]; b_val <- stats::coef(m)["b"]
      dma_eq <- calc_dma_especifico(dados_st, "Eq5", a_val, b_val)
      resultados_st <- rbind(resultados_st, data.frame(
        amostra_id = id, Equacao = "Eq5", R2 = calc_r2(F_obs, stats::predict(m)),
        DMA_calculado = dma_eq, Param_a = a_val, Param_b = b_val
      ))
    }, error = function(e) {
      resultados_st <<- rbind(resultados_st, data.frame(amostra_id = id, Equacao = "Eq5", R2 = NA, DMA_calculado = NA, Param_a = NA, Param_b = NA))
    })

    lista_detalhes_modelos[[i]] <- resultados_st

    # Selecionar a melhor (pelo R2) para o resumo final
    # Remove NAs para a comparação
    res_validos <- resultados_st[!is.na(resultados_st$R2), ]
    if(nrow(res_validos) > 0) {
      melhor <- res_validos[which.max(res_validos$R2), ]
      lista_resumo_dma[[i]] <- data.frame(
        amostra_id = id, DMA = melhor$DMA_calculado,
        Melhor_Equacao = melhor$Equacao, R2_Melhor = melhor$R2
      )
    } else {
      lista_resumo_dma[[i]] <- data.frame(amostra_id = id, DMA = NA, Melhor_Equacao = NA, R2_Melhor = NA)
    }
  }

  # Consolidar tabelas
  tab_dma <- do.call(rbind, lista_resumo_dma)
  tab_detalhes <- do.call(rbind, lista_detalhes_modelos)

  # Renomear colunas para bater com o input original do usuário
  colnames(tab_dma)[1] <- col_amostra
  colnames(tab_detalhes)[1] <- col_amostra

  # Cruzando com DMP e DMG
  resumo_final <- dplyr::left_join(indices_trad, tab_dma, by = col_amostra)

  # Adicionando a coluna de validação da recomendação dos autores
  resumo_final$Uso_Recomendado <- ifelse(resumo_final$R2_Melhor >= 0.99, "Sim", "Não (R2 < 0.99)")

  # Verificação da regra de van Lier & Albuquerque (1997) para disparar o aviso
  if (any(resumo_final$R2_Melhor < 0.99, na.rm = TRUE)) {
    warning("Atenção: Pelo menos uma amostra apresentou R2 máximo inferior a 0,99. ",
            "van Lier & Albuquerque (1997) recomendam o uso do DMA apenas quando R2 >= 0,99. ",
            "Verifique a coluna 'Uso_Recomendado' no resumo gerado.")
  }

  return(list(
    resumo = resumo_final,
    detalhes_equacoes = tab_detalhes
  ))
}


#' Gráfico de diagnóstico comparativo para uma amostra
#'
#' Gera uma figura com três gráficos (Eq3, Eq4 e Eq5) para uma amostra específica,
#' exibindo o R2, os parâmetros de ajuste e o DMA calculado em cada um.
#'
#' @note Segundo van Lier & Albuquerque (1997), o método proposto para o cálculo
#' do DMA deve ser aplicado sem ressalvas apenas quando o coeficiente de
#' determinação (R2) do ajuste da equação for maior ou igual a 0,99.
#'
#' @param dados Data frame bruto no formato longo.
#' @param amostra_id_alvo O identificador da amostra que deseja analisar.
#' @param col_amostra Nome da coluna de amostra.
#' @param col_diametro Nome da coluna de diâmetro (mm).
#' @param col_massa Nome da coluna de massa (g).
#'
#' @export
#'
plot_diagnostico_amostra <- function(dados, amostra_id_alvo, col_amostra = "amostra_id",
                                     col_diametro = "diametro_mm", col_massa = "massa_g") {

  # 1. Filtrar e preparar dados da amostra
  dados_st <- dados[dados[[col_amostra]] == amostra_id_alvo, ]
  prep <- prep_agregados(dados_st, col_amostra, col_diametro, col_massa)
  df_proc <- prep$dados_processados

  D <- df_proc[[col_diametro]]
  F_obs <- df_proc$F_acumulada
  soma_F <- df_proc$F_acumulada + df_proc$F_acumulada_inf

  # Funções auxiliares internas
  calc_dma_val <- function(eq, a, b) {
    if (eq == "Eq3") D_aj <- (b * soma_F) / (2 - (a * soma_F))
    else if (eq == "Eq4") D_aj <- ((b * soma_F) / (2 - (a * soma_F)))^2
    else D_aj <- (soma_F / (2 * a))^(1 / b)
    D_aj[!is.finite(D_aj)] <- 0
    return(sum(df_proc$fracao_massa * D_aj, na.rm = TRUE))
  }

  calc_r2 <- function(obs, pred) {
    sq_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    sq_res <- sum((obs - pred)^2, na.rm = TRUE)
    return(1 - (sq_res / sq_tot))
  }

  # 2. Processar as 3 equações para gerar dados de curva e labels
  plot_data_list <- list()
  curva_data_list <- list()

  equacoes <- c("Eq3", "Eq4", "Eq5")
  formulas <- list(
    Eq3 = function(D, a, b) 1 / (a + (b / D)),
    Eq4 = function(D, a, b) 1 / (a + (b / sqrt(D))),
    Eq5 = function(D, a, b) a * (D^b)
  )

  for (eq in equacoes) {
    a_val <- NA; b_val <- NA; dma_val <- NA; r2_val <- NA

    tryCatch({
      if(eq == "Eq3") m <- stats::nls(F_obs ~ 1 / (a + (b / D)), start = list(a = 1, b = 1))
      else if(eq == "Eq4") m <- stats::nls(F_obs ~ 1 / (a + (b / sqrt(D))), start = list(a = 1, b = 1))
      else m <- stats::nls(F_obs ~ a * (D^b), start = list(a = 1, b = 1))

      a_val <- stats::coef(m)["a"]; b_val <- stats::coef(m)["b"]
      dma_val <- calc_dma_val(eq, a_val, b_val)
      r2_val <- calc_r2(F_obs, stats::predict(m))

      # Gerar curva teórica
      d_seq <- seq(0.01, max(D), length.out = 100)
      f_seq <- formulas[[eq]](d_seq, a_val, b_val)
      curva_data_list[[eq]] <- data.frame(diametro_mm = d_seq, F_teorica = f_seq, Equacao = eq)

    }, error = function(e) { })

    # Criar label informativo para o facet (Agora com R2)
    label_info <- paste0(eq, "\n",
                         "R² = ", ifelse(is.na(r2_val), "NA", round(r2_val, 4)), "\n",
                         "a = ", ifelse(is.na(a_val), "NA", round(a_val, 3)), "\n",
                         "b = ", ifelse(is.na(b_val), "NA", round(b_val, 3)), "\n",
                         "DMA = ", ifelse(is.na(dma_val), "NA", round(dma_val, 3)), " mm")

    # Dados dos pontos observados com o label da faceta
    plot_data_list[[eq]] <- data.frame(
      diametro_mm = D,
      F_acumulada = F_obs,
      Equacao_Label = label_info,
      DMA_ref = dma_val
    )

    if (!is.null(curva_data_list[[eq]])) {
      curva_data_list[[eq]]$Equacao_Label <- label_info
    }
  }

  df_pontos <- do.call(rbind, plot_data_list)
  df_curvas <- do.call(rbind, curva_data_list)

  # 3. Gerar o gráfico
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_curvas, ggplot2::aes(x = diametro_mm, y = F_teorica), color = "blue") +
    ggplot2::geom_point(data = df_pontos, ggplot2::aes(x = diametro_mm, y = F_acumulada)) +
    ggplot2::geom_vline(data = df_pontos, ggplot2::aes(xintercept = DMA_ref), color = "red", linetype = "dashed") +
    ggplot2::facet_wrap(~Equacao_Label) +
    ggplot2::labs(title = paste("Diagnóstico de Modelos - Amostra:", amostra_id_alvo),
                  x = "Diâmetro (mm)", y = "Fração Acumulada (g/g)",
                  caption = "Nota: van Lier & Albuquerque (1997) recomendam o uso do DMA apenas se R² >= 0,99.") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text = ggplot2::element_text(face = "bold", size = 10)
    )

  return(p)
}
