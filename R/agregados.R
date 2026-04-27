#' Prepara dados de peneiramento e calcula DMP e DMG
#'
#' @export
prep_agregados <- function(df, col_amostra = "amostra_id", col_diametro = "diametro_mm", col_massa = "massa_g") {

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
#' @export
calc_dma <- function(df_processado, col_amostra = "amostra_id", col_diametro = "diametro_mm") {

  calcular_r2 <- function(observado, predito) {
    sq_tot <- sum((observado - mean(observado, na.rm = TRUE))^2, na.rm = TRUE)
    sq_res <- sum((observado - predito)^2, na.rm = TRUE)
    return(1 - (sq_res / sq_tot))
  }

  amostras_lista <- split(df_processado, df_processado[[col_amostra]])

  resultados_dma <- lapply(amostras_lista, function(dados_amostra) {
    nome_amostra <- unique(dados_amostra[[col_amostra]])

    D_vetor <- dados_amostra[[col_diametro]]
    F_vetor <- dados_amostra$F_acumulada

    modelos <- list()

    tryCatch({
      mod1 <- stats::nls(F_vetor ~ 1 / (a + (b / D_vetor)), start = list(a = 1, b = 1))
      modelos[["Eq3"]] <- list(nome = "Eq3", a = stats::coef(mod1)["a"], b = stats::coef(mod1)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod1)))
    }, error = function(e) NULL)

    tryCatch({
      mod2 <- stats::nls(F_vetor ~ 1 / (a + (b / sqrt(D_vetor))), start = list(a = 1, b = 1))
      modelos[["Eq4"]] <- list(nome = "Eq4", a = stats::coef(mod2)["a"], b = stats::coef(mod2)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod2)))
    }, error = function(e) NULL)

    tryCatch({
      mod3 <- stats::nls(F_vetor ~ a * (D_vetor^b), start = list(a = 1, b = 1))
      modelos[["Eq5"]] <- list(nome = "Eq5", a = stats::coef(mod3)["a"], b = stats::coef(mod3)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod3)))
    }, error = function(e) NULL)

    tryCatch({
      mod4 <- stats::nls(F_vetor ~ 1 - exp(-(D_vetor/a)^b), start = list(a = max(D_vetor)/2, b = 1))
      modelos[["Weibull"]] <- list(nome = "Weibull", a = stats::coef(mod4)["a"], b = stats::coef(mod4)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod4)))
    }, error = function(e) NULL)

    tryCatch({
      mod5 <- stats::nls(F_vetor ~ 1 / (1 + exp(-a * (D_vetor - b))), start = list(a = 1, b = mean(D_vetor)))
      modelos[["Logistic"]] <- list(nome = "Logistic", a = stats::coef(mod5)["a"], b = stats::coef(mod5)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod5)))
    }, error = function(e) NULL)

    tryCatch({
      mod6 <- stats::nls(F_vetor ~ stats::plnorm(D_vetor, a, b), start = list(a = 0, b = 1))
      modelos[["LogNormal"]] <- list(nome = "LogNormal", a = stats::coef(mod6)["a"], b = stats::coef(mod6)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod6)))
    }, error = function(e) NULL)

    tryCatch({
      mod7 <- stats::nls(F_vetor ~ exp(-exp(-a * (D_vetor - b))), start = list(a = 1, b = mean(D_vetor)))
      modelos[["Gompertz"]] <- list(nome = "Gompertz", a = stats::coef(mod7)["a"], b = stats::coef(mod7)["b"], r2 = calcular_r2(F_vetor, stats::predict(mod7)))
    }, error = function(e) NULL)

    modelos <- Filter(Negate(is.null), modelos)

    if (length(modelos) == 0) {
      return(data.frame(amostra_id = nome_amostra, DMA = NA, Modelo = NA, R2 = NA, a = NA, b = NA))
    }

    melhor <- modelos[[which.max(sapply(modelos, function(x) x$r2))]]

    soma_F <- dados_amostra$F_acumulada + dados_amostra$F_acumulada_inf
    F_m <- soma_F / 2

    if (melhor$nome == "Eq3") {
      D_ajustado <- (melhor$b * soma_F) / (2 - (melhor$a * soma_F))
    } else if (melhor$nome == "Eq4") {
      D_ajustado <- ((melhor$b * soma_F) / (2 - (melhor$a * soma_F)))^2
    } else if (melhor$nome == "Eq5") {
      D_ajustado <- (soma_F / (2 * melhor$a))^(1 / melhor$b)
    } else if (melhor$nome == "Weibull") {
      D_ajustado <- melhor$a * (-log(1 - pmin(F_m, 0.9999)))^(1 / melhor$b)
    } else if (melhor$nome == "Logistic") {
      F_m_adj <- pmax(pmin(F_m, 0.9999), 0.0001)
      D_ajustado <- melhor$b - (1 / melhor$a) * log((1 / F_m_adj) - 1)
    } else if (melhor$nome == "LogNormal") {
      F_m_adj <- pmax(pmin(F_m, 0.9999), 0.0001)
      D_ajustado <- stats::qlnorm(F_m_adj, melhor$a, melhor$b)
    } else if (melhor$nome == "Gompertz") {
      F_m_adj <- pmax(pmin(F_m, 0.9999), 0.0001)
      D_ajustado <- melhor$b - (1 / melhor$a) * log(-log(F_m_adj))
    }

    D_ajustado[!is.finite(D_ajustado)] <- 0
    D_ajustado[D_ajustado < 0] <- 0
    D_ponderado <- dados_amostra$fracao_massa * D_ajustado
    DMA_final <- sum(D_ponderado, na.rm = TRUE)

    df_resumo <- data.frame(
      amostra_id = nome_amostra,
      DMA = DMA_final,
      Modelo_Vencedor = melhor$nome,
      R2 = melhor$r2,
      Param_a = melhor$a,
      Param_b = melhor$b
    )

    colnames(df_resumo)[1] <- col_amostra
    return(df_resumo)
  })

  tabela_final <- do.call(rbind, resultados_dma)
  rownames(tabela_final) <- NULL

  if (any(tabela_final$R2 < 0.99, na.rm = TRUE)) {
    warning("Atenção: Pelo menos uma amostra apresentou R2 máximo inferior a 0,99. ",
            "van Lier & Albuquerque (1997) recomendam o uso do DMA apenas quando R2 >= 0,99.")
  }

  return(tabela_final)
}


#' Plota as curvas de retenção e o ajuste do DMA
#'
#' @export
plot_dma <- function(df_processado, df_dma, col_amostra = "amostra_id", col_diametro = "diametro_mm") {

  gerar_curva <- function(amostra, mod, param_a, param_b, d_max) {
    d_seq <- seq(0.01, d_max, length.out = 100)

    if (mod == "Eq3") {
      f_seq <- 1 / (param_a + (param_b / d_seq))
    } else if (mod == "Eq4") {
      f_seq <- 1 / (param_a + (param_b / sqrt(d_seq)))
    } else if (mod == "Eq5") {
      f_seq <- param_a * (d_seq^param_b)
    } else if (mod == "Weibull") {
      f_seq <- 1 - exp(-(d_seq/param_a)^param_b)
    } else if (mod == "Logistic") {
      f_seq <- 1 / (1 + exp(-param_a * (d_seq - param_b)))
    } else if (mod == "LogNormal") {
      f_seq <- stats::plnorm(d_seq, param_a, param_b)
    } else if (mod == "Gompertz") {
      f_seq <- exp(-exp(-param_a * (d_seq - param_b)))
    }

    df <- data.frame(amostra_id = amostra, diametro_mm = d_seq, F_teorica = f_seq)
    colnames(df)[1] <- col_amostra
    colnames(df)[2] <- col_diametro
    return(df)
  }

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

  df_curvas <- do.call(rbind, lista_curvas)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_curvas, ggplot2::aes(x = .data[[col_diametro]], y = F_teorica), color = "blue", linewidth = 0.8) +
    ggplot2::geom_point(data = df_processado, ggplot2::aes(x = .data[[col_diametro]], y = F_acumulada), size = 2) +
    ggplot2::geom_vline(data = df_dma, ggplot2::aes(xintercept = DMA), color = "darkred", linetype = "dashed") +
    ggplot2::facet_wrap(stats::as.formula(paste("~", col_amostra))) +
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
      strip.text = ggplot2::element_text(face = "bold")
    )

  return(p)
}


#' Executa análise completa e detalhada de agregação
#'
#' @export
analise_completa_dma <- function(dados, col_amostra = "amostra_id", col_diametro = "diametro_mm", col_massa = "massa_g") {

  prep <- prep_agregados(dados, col_amostra, col_diametro, col_massa)
  df_proc <- prep$dados_processados
  indices_trad <- prep$indices

  calc_r2 <- function(obs, pred) {
    sq_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    sq_res <- sum((obs - pred)^2, na.rm = TRUE)
    return(1 - (sq_res / sq_tot))
  }

  calc_dma_especifico <- function(dados_st, eq_nome, a, b) {
    soma_F <- dados_st$F_acumulada + dados_st$F_acumulada_inf
    F_m <- soma_F / 2

    if (eq_nome == "Eq3") {
      D_aj <- (b * soma_F) / (2 - (a * soma_F))
    } else if (eq_nome == "Eq4") {
      D_aj <- ((b * soma_F) / (2 - (a * soma_F)))^2
    } else if (eq_nome == "Eq5") {
      D_aj <- (soma_F / (2 * a))^(1 / b)
    } else if (eq_nome == "Weibull") {
      D_aj <- a * (-log(1 - pmin(F_m, 0.9999)))^(1 / b)
    } else if (eq_nome == "Logistic") {
      F_m_adj <- pmax(pmin(F_m, 0.9999), 0.0001)
      D_aj <- b - (1 / a) * log((1 / F_m_adj) - 1)
    } else if (eq_nome == "LogNormal") {
      F_m_adj <- pmax(pmin(F_m, 0.9999), 0.0001)
      D_aj <- stats::qlnorm(F_m_adj, a, b)
    } else if (eq_nome == "Gompertz") {
      F_m_adj <- pmax(pmin(F_m, 0.9999), 0.0001)
      D_aj <- b - (1 / a) * log(-log(F_m_adj))
    }

    D_aj[!is.finite(D_aj)] <- 0
    D_aj[D_aj < 0] <- 0
    return(sum(dados_st$fracao_massa * D_aj, na.rm = TRUE))
  }

  amostras_lista <- split(df_proc, df_proc[[col_amostra]])
  lista_resumo_dma <- list()
  lista_detalhes_modelos <- list()

  tentar_ajuste <- function(formula_nls, start_list, eq_nome, D_vetor, F_obs, dados_st, id) {
    tryCatch({
      m <- stats::nls(formula_nls, start = start_list)
      a_val <- stats::coef(m)["a"]; b_val <- stats::coef(m)["b"]
      dma_eq <- calc_dma_especifico(dados_st, eq_nome, a_val, b_val)
      data.frame(
        amostra_id = id, Equacao = eq_nome, R2 = calc_r2(F_obs, stats::predict(m)),
        DMA_calculado = dma_eq, Param_a = a_val, Param_b = b_val
      )
    }, error = function(e) {
      data.frame(amostra_id = id, Equacao = eq_nome, R2 = NA, DMA_calculado = NA, Param_a = NA, Param_b = NA)
    })
  }

  for (i in seq_along(amostras_lista)) {
    dados_st <- amostras_lista[[i]]
    id <- unique(dados_st[[col_amostra]])
    D <- dados_st[[col_diametro]]
    F_obs <- dados_st$F_acumulada

    res_eq3 <- tentar_ajuste(F_obs ~ 1 / (a + (b / D)), list(a = 1, b = 1), "Eq3", D, F_obs, dados_st, id)
    res_eq4 <- tentar_ajuste(F_obs ~ 1 / (a + (b / sqrt(D))), list(a = 1, b = 1), "Eq4", D, F_obs, dados_st, id)
    res_eq5 <- tentar_ajuste(F_obs ~ a * (D^b), list(a = 1, b = 1), "Eq5", D, F_obs, dados_st, id)
    res_weibull <- tentar_ajuste(F_obs ~ 1 - exp(-(D/a)^b), list(a = max(D)/2, b = 1), "Weibull", D, F_obs, dados_st, id)
    res_logistic <- tentar_ajuste(F_obs ~ 1 / (1 + exp(-a * (D - b))), list(a = 1, b = mean(D)), "Logistic", D, F_obs, dados_st, id)
    res_lognormal <- tentar_ajuste(F_obs ~ stats::plnorm(D, a, b), list(a = 0, b = 1), "LogNormal", D, F_obs, dados_st, id)
    res_gompertz <- tentar_ajuste(F_obs ~ exp(-exp(-a * (D - b))), list(a = 1, b = mean(D)), "Gompertz", D, F_obs, dados_st, id)

    resultados_st <- rbind(res_eq3, res_eq4, res_eq5, res_weibull, res_logistic, res_lognormal, res_gompertz)
    lista_detalhes_modelos[[i]] <- resultados_st

    res_validos <- resultados_st[!is.na(resultados_st$R2), ]
    if(nrow(res_validos) > 0) {
      melhor <- res_validos[which.max(res_validos$R2), ]

      # --- A CORREÇÃO ESTÁ AQUI NESTE BLOCO ---
      lista_resumo_dma[[i]] <- data.frame(
        amostra_id = id,
        DMA = melhor$DMA_calculado,
        Modelo_Vencedor = melhor$Equacao, # Nome corrigido para combinar com o plot
        R2_Melhor = melhor$R2,
        Param_a = melhor$Param_a,         # Adicionado para o plot
        Param_b = melhor$Param_b          # Adicionado para o plot
      )
    } else {
      lista_resumo_dma[[i]] <- data.frame(
        amostra_id = id, DMA = NA, Modelo_Vencedor = NA,
        R2_Melhor = NA, Param_a = NA, Param_b = NA
      )
    }
  }

  tab_dma <- do.call(rbind, lista_resumo_dma)
  tab_detalhes <- do.call(rbind, lista_detalhes_modelos)
  colnames(tab_dma)[1] <- col_amostra
  colnames(tab_detalhes)[1] <- col_amostra

  resumo_final <- dplyr::left_join(indices_trad, tab_dma, by = col_amostra)
  resumo_final$Uso_Recomendado <- ifelse(resumo_final$R2_Melhor >= 0.99, "Sim", "Não (R2 < 0.99)")

  if (any(resumo_final$R2_Melhor < 0.99, na.rm = TRUE)) {
    warning("Atenção: Pelo menos uma amostra apresentou R2 máximo inferior a 0,99. ",
            "van Lier & Albuquerque (1997) recomendam o uso do DMA apenas quando R2 >= 0,99.")
  }

  return(list(resumo = resumo_final, detalhes_equacoes = tab_detalhes))
}


#' Gráfico de diagnóstico comparativo para uma amostra
#'
#' @export
plot_diagnostico_amostra <- function(dados, amostra_id_alvo, col_amostra = "amostra_id",
                                     col_diametro = "diametro_mm", col_massa = "massa_g") {

  dados_st <- dados[dados[[col_amostra]] == amostra_id_alvo, ]
  prep <- prep_agregados(dados_st, col_amostra, col_diametro, col_massa)
  df_proc <- prep$dados_processados

  D <- df_proc[[col_diametro]]
  F_obs <- df_proc$F_acumulada
  soma_F <- df_proc$F_acumulada + df_proc$F_acumulada_inf
  F_m <- soma_F / 2

  calc_dma_val <- function(eq, a, b) {
    if (eq == "Eq3") D_aj <- (b * soma_F) / (2 - (a * soma_F))
    else if (eq == "Eq4") D_aj <- ((b * soma_F) / (2 - (a * soma_F)))^2
    else if (eq == "Eq5") D_aj <- (soma_F / (2 * a))^(1 / b)
    else if (eq == "Weibull") D_aj <- a * (-log(1 - pmin(F_m, 0.9999)))^(1 / b)
    else if (eq == "Logistic") D_aj <- b - (1 / a) * log((1 / pmax(pmin(F_m, 0.9999), 0.0001)) - 1)
    else if (eq == "LogNormal") D_aj <- stats::qlnorm(pmax(pmin(F_m, 0.9999), 0.0001), a, b)
    else if (eq == "Gompertz") D_aj <- b - (1 / a) * log(-log(pmax(pmin(F_m, 0.9999), 0.0001)))

    D_aj[!is.finite(D_aj)] <- 0
    D_aj[D_aj < 0] <- 0
    return(sum(df_proc$fracao_massa * D_aj, na.rm = TRUE))
  }

  calc_r2 <- function(obs, pred) {
    sq_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    sq_res <- sum((obs - pred)^2, na.rm = TRUE)
    return(1 - (sq_res / sq_tot))
  }

  plot_data_list <- list()
  curva_data_list <- list()

  equacoes <- c("Eq3", "Eq4", "Eq5", "Weibull", "Logistic", "LogNormal", "Gompertz")
  formulas <- list(
    Eq3 = function(D, a, b) 1 / (a + (b / D)),
    Eq4 = function(D, a, b) 1 / (a + (b / sqrt(D))),
    Eq5 = function(D, a, b) a * (D^b),
    Weibull = function(D, a, b) 1 - exp(-(D/a)^b),
    Logistic = function(D, a, b) 1 / (1 + exp(-a * (D - b))),
    LogNormal = function(D, a, b) stats::plnorm(D, a, b),
    Gompertz = function(D, a, b) exp(-exp(-a * (D - b)))
  )

  for (eq in equacoes) {
    a_val <- NA; b_val <- NA; dma_val <- NA; r2_val <- NA

    tryCatch({
      if(eq == "Eq3") m <- stats::nls(F_obs ~ 1 / (a + (b / D)), start = list(a = 1, b = 1))
      else if(eq == "Eq4") m <- stats::nls(F_obs ~ 1 / (a + (b / sqrt(D))), start = list(a = 1, b = 1))
      else if(eq == "Eq5") m <- stats::nls(F_obs ~ a * (D^b), start = list(a = 1, b = 1))
      else if(eq == "Weibull") m <- stats::nls(F_obs ~ 1 - exp(-(D/a)^b), start = list(a = max(D)/2, b = 1))
      else if(eq == "Logistic") m <- stats::nls(F_obs ~ 1 / (1 + exp(-a * (D - b))), start = list(a = 1, b = mean(D)))
      else if(eq == "LogNormal") m <- stats::nls(F_obs ~ stats::plnorm(D, a, b), start = list(a = 0, b = 1))
      else if(eq == "Gompertz") m <- stats::nls(F_obs ~ exp(-exp(-a * (D - b))), start = list(a = 1, b = mean(D)))

      a_val <- stats::coef(m)["a"]; b_val <- stats::coef(m)["b"]
      dma_val <- calc_dma_val(eq, a_val, b_val)
      r2_val <- calc_r2(F_obs, stats::predict(m))

      d_seq <- seq(0.01, max(D), length.out = 100)
      f_seq <- formulas[[eq]](d_seq, a_val, b_val)
      curva_data_list[[eq]] <- data.frame(diametro_mm = d_seq, F_teorica = f_seq, Equacao = eq)

    }, error = function(e) { })

    label_info <- paste0(eq, "\n",
                         "R² = ", ifelse(is.na(r2_val), "NA", round(r2_val, 4)), "\n",
                         "a = ", ifelse(is.na(a_val), "NA", round(a_val, 3)), "\n",
                         "b = ", ifelse(is.na(b_val), "NA", round(b_val, 3)), "\n",
                         "DMA = ", ifelse(is.na(dma_val), "NA", round(dma_val, 3)), " mm")

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

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_curvas, ggplot2::aes(x = diametro_mm, y = F_teorica), color = "blue") +
    ggplot2::geom_point(data = df_pontos, ggplot2::aes(x = diametro_mm, y = F_acumulada)) +
    ggplot2::geom_vline(data = df_pontos, ggplot2::aes(xintercept = DMA_ref), color = "red", linetype = "dashed") +
    ggplot2::facet_wrap(~Equacao_Label, ncol = 3) +
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

#' Plota a distribuição de tamanho de agregados
#'
#' Cria um gráfico de barras ilustrando a proporção de massa retida
#' em cada classe de diâmetro para cada amostra.
#'
#' @param df_processado Data frame. O objeto \code{dados_processados} retornado pela função \code{prep_agregados()}.
#' @param col_amostra String. Nome da coluna de identificação da amostra.
#' @param col_diametro String. Nome da coluna de diâmetro limite das classes (mm).
#'
#' @return Um objeto ggplot com o painel de distribuição por amostra.
#' @export
plot_distribuicao_agregados <- function(df_processado, col_amostra = "amostra_id", col_diametro = "diametro_mm") {

  # Criar uma cópia isolada para plotagem
  df_plot <- df_processado

  # Converter fração para porcentagem para o eixo Y
  df_plot$Porcentagem <- df_plot$fracao_massa * 100

  # Ordenar o eixo X (peneiras) de forma decrescente (do maior diâmetro para o menor)
  niveis_ordem <- sort(unique(df_plot[[col_diametro]]), decreasing = TRUE)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(
    x = factor(.data[[col_diametro]], levels = niveis_ordem),
    y = Porcentagem
  )) +
    ggplot2::geom_col(fill = "#2c3e50", color = "black", alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = round(Porcentagem, 1)),
                       vjust = -0.5, size = 3, color = "black") +
    ggplot2::facet_wrap(stats::as.formula(paste("~", col_amostra))) +
    ggplot2::labs(
      title = "Distribuição de Agregados por Classe de Tamanho",
      x = "Diâmetro Limite da Classe (mm)",
      y = "Massa Retida (%)",
      caption = "Valores sobre as barras representam a porcentagem de massa."
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # Adiciona uma margem superior para o texto não cortar
  p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))

  return(p)
}



#' Exporta os resultados da análise de agregação para Excel
#'
#' Esta função recebe o objeto de lista gerado pela função
#' \code{analise_completa_dma()} e salva as tabelas de resumo e
#' detalhes em abas separadas de um arquivo .xlsx.
#'
#' @param lista_analise Lista. O objeto retornado por \code{analise_completa_dma()}.
#' @param caminho_arquivo String. O nome ou caminho completo do arquivo a ser criado.
#'
#' @return A função não retorna um objeto no R, mas salva o arquivo no disco.
#' @export
#'
exportar_analise_xlsx <- function(lista_analise, caminho_arquivo = "resultados_dma.xlsx") {

  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop(
      "O pacote 'writexl' é necessário para esta função. ",
      "Por favor, instale-o com install.packages('writexl')"
    )
  }

  if (!is.list(lista_analise) || !all(c("resumo", "detalhes_equacoes") %in% names(lista_analise))) {
    stop("A entrada deve ser a lista retornada pela função analise_completa_dma().")
  }

  dados_para_exportar <- list(
    "Resumo_Geral" = lista_analise$resumo,
    "Detalhes_dos_Modelos" = lista_analise$detalhes_equacoes
  )

  writexl::write_xlsx(dados_para_exportar, path = caminho_arquivo)

  message(paste("Arquivo exportado com sucesso para:", caminho_arquivo))
}
