#> heckitmfx()                                               # 限界効果を推定するジェネリック関数（モデルのクラスでメソッド分派）
#> ├─ heckitmfx.selection()                                  # sampleSelection::heckit 用メソッド：係数抽出・設計行列Z作成→計算
#> │  ├─ coef(), model.matrix(), stringr::str_remove()       # 係数/設計行列の取得と前処理
#> │  └─ heckitmfx_compute()                                 # 与えた係数とZでME（無条件/条件/選択）を一括計算
#> │     ├─ is_dummy()                                       # Zの各列が {0,1} のダミーか判定
#> │     ├─ terms_for_dummy()                                # ダミーを0→1に切替えた時の間接効果（IMR差・選択確率の対数差）
#> │     └─ （dnorm, pnorm 等の基礎関数・dplyr操作）         # 逆ミルズ比λ、δの計算とテーブル整形
#> └─ heckitmfx.heckitrob()                                  # ssmrob::heckitrob 用メソッド：係数抽出・Z作成→計算
#> ├─ coef(), vcov(), stringr::str_remove(), model.matrix()  # 係数/設計行列の取得と前処理（robust版）
#> └─ heckitmfx_compute()                                    # 上と同じ計算エンジンを利用
#> ├─ is_dummy()
#> ├─ terms_for_dummy()
#> └─ （dnorm, pnorm 等・dplyr操作）
#> 
#> heckitmfx_delta()                                         # デルタ法で標準誤差・Z値・P値を付与
#> ├─ heckitmfx()                                            # まず推定量（unconditional/conditional/selection のいずれか）を取得
#> ├─ numDeriv::jacobian()                                   # 係数に対するヤコビアン（数値微分）を計算
#> ├─ vcov() / checkmate::testClass()                        # パラメータ共分散行列の構築（heckitrobは段結合，heckitはNA→0）
#> └─ 行列演算（grad %*% V %*% t(grad) → se）                 # デルタ法で標準誤差→統計量・P値を算出
#> 
#> is_dummy()                                                # ベクトルが {0,1} だけかを判定（checkmate::test_set_equal を使用）
#> terms_for_dummy()                                         # （Vectorize）ダミー0/1でのIMR差×βλ と logΦ差を返す（行列[2×#dummy]）
#> heckitmfx_compute()                                       # 係数とZから ME の各項を計算し、整形して返す中核関数

# メインの関数から呼び出されて、限界効果を推定する関数 ------------------------------------------------------------
heckitmfx_compute <- function(
    gamma, beta, betalambda, Z, ...
    # type = c("unconditional", "conditional", "selection")
){
  # 就業決定関数の予測値と逆ミルズ比の計算 ---------------------------------
  alpha <- Z %*% gamma
  lambda <- dnorm(alpha) / pnorm(alpha)
  delta <- lambda * (lambda + alpha)
  
  # 限界効果の各項の計算 ---------------------
  selection <- gamma[-1] * mean(lambda)
  ei_2 <- gamma[-1] * betalambda * mean(delta)
  
  # ダミー変数用の処理 --------------
  dummy_vars <- apply(Z[, -1], FUN = is_dummy, MARGIN = 2)
  
  if(sum(dummy_vars) >= 1){
    terms_for_dummy_vars <- names(dummy_vars[dummy_vars]) |>
      terms_for_dummy(Z, gamma, betalambda)
    
    # ダミー変数についての推定値に置換
    ei_2[dummy_vars] <- terms_for_dummy_vars[1, ]
    selection[dummy_vars] <- terms_for_dummy_vars[2, ]
  }
  
  # 回帰係数と限界効果の各項をデータフレームにまとめます。
  est <- merge(
    data.frame(beta = beta[-1]), 
    data.frame(gamma = gamma[-1]), 
    by = "row.names", all = TRUE
  ) |> 
    rename(term = "Row.names") |> 
    full_join(tibble::tibble(ei_2, term = names(ei_2)), by = join_by(term)) |> 
    full_join(tibble::tibble(selection, term = names(selection)), by = join_by(term))
  
  est[is.na(est)] <- 0
  
  # 限界効果の推定
  res <- est |> 
    mutate(
      term = as.character(term),
      conditional = beta - ei_2,
      unconditional = conditional + selection
    ) |> 
    select(term, unconditional, conditional, selection, beta, gamma)
  # 結果の出力 ---------------------------------
  
  return(res)
}


# ダミー変数を判定する関数
is_dummy <- function(x, unique_value = c(0, 1), na.rm = FALSE) {
  if(na.rm) x <- na.omit(x)
  res <- checkmate::test_set_equal(x = x, y = unique_value)
  return(res)
}

# ダミー変数に関する間接効果の計算
terms_for_dummy <- Vectorize(function(var_name, Z, gamma, betalambda){
  fun_inv_mills <- \(x) dnorm(x) / pnorm(x)
  
  z_1 <- z_0 <- colMeans(Z)
  z_1[var_name] <- 1
  z_0[var_name] <- 0
  
  delta_lambda <- betalambda * as.vector(fun_inv_mills(z_1 %*% gamma) - fun_inv_mills(z_0 %*% gamma))
  delta_log_phi <- as.vector(pnorm(z_1 %*% gamma, log.p = TRUE) - pnorm(z_0 %*% gamma, log.p = TRUE))
  
  res <- rbind(delta_lambda, delta_log_phi) 
  return(res)
}, vectorize.args = "var_name")


# モデルオブジェクトから限界効果を推定するメインの関数 ------------------------------------------------------------
# ジェネリック関数の定義
heckitmfx <- function(model, newdata = NULL, .params = NULL,...) UseMethod("heckitmfx")

# sampleSelection::heckit() 用のメソッド -----------------
heckitmfx.selection <- function(
    model, newdata = NULL, .params = NULL, ...
){
  # 回帰係数の抽出 ---------------------------------
  if(is.null(.params)){
    coef_all <- coef(model) # 全ての回帰係数を取得  
  }else{
    # デルタ法の実装に使うための、回帰係数を任意の値で上書きする機能
    coef_all <- .params # 引数 .params が指定されていた場合は .params で上書きする
  }
  
  # 全ての説明変数の名前を抽出
  all_vars <- names(coef_all) |> unique() %>% 
    subset(!(. %in% c("(Intercept)", "invMillsRatio", "sigma", "rho")))
  
  # 賃金関数の回帰係数が始まる位置を特定する
  start_beta <- which(names(coef_all) == "(Intercept)")[2]
  
  gamma <- coef_all[1:(start_beta - 1)] # 就業決定関数の回帰係数
  # return(res_gamma)
  
  # 賃金関数の回帰係数
  beta <- coef_all[start_beta:length(coef_all)] %>%
    subset(!(names(.) %in% c("invMillsRatio", "sigma", "rho")))
  
  betalambda <- coef_all %>% subset((names(.) %in% c("sigma", "rho"))) |> prod()
  
  # 就業決定関数の計画行列の抽出 ---------------------------------
  if(model$method == "2step"){
    Z <- model.matrix(model$probit, data = newdata)
    colnames(Z) <- str_remove(colnames(Z), "XS")
  }else{
    if(model$method == "ml"){
      Z <- model.matrix(model$termsS, data = newdata)
    }
  }
  # 限界効果の推定
  res <- heckitmfx_compute(
    gamma = gamma,
    beta = beta,
    betalambda = betalambda,
    Z = Z
  )
  return(res)
}

# ssmrob::heckitrob() 用のメソッド -----------------
heckitmfx.heckitrob <- function(
    model, newdata = NULL, .params = NULL, ...
){
  # 回帰係数の抽出 ---------------------------------
  if(is.null(.params)){
    # 就業決定関数の回帰係数
    gamma <- coef(model$stage1)
    # 賃金関数の回帰係数
    beta <- coef(model$stage2) %>% subset(names(.) != "IMR1")
    betalambda <- coef(model$stage2)["IMR1"]
  }else{
    # デルタ法の実装に使うための、回帰係数を任意の値で上書きする機能
    coef_all <- .params # 引数 .params が指定されていた場合は .params で上書きする
    names(coef_all) <- str_remove(names(coef_all), "S.|O.")
    
    # return(coef_all)
    # 全ての説明変数の名前を抽出
    all_vars <- names(coef_all) |> unique() %>% 
      subset(!(. %in% c("(Intercept)", "IMR1")))
    
    # 賃金関数の回帰係数が始まる位置を特定する
    start_beta <- which(names(coef_all) == "(Intercept)")[2]
    
    gamma <- coef_all[1:(start_beta - 1)] # 就業決定関数の回帰係数
    
    # 賃金関数の回帰係数
    beta <- coef_all[start_beta:length(coef_all)] %>% subset(names(.) != "IMR1")
    betalambda <- coef_all["IMR1"]
  }
  
  # 就業決定関数の計画行列の抽出 ---------------------------------
  Z <- model.matrix(model$stage1, data = newdata)
  # 限界効果の推定
  res <- heckitmfx_compute(
    gamma = gamma,
    beta = beta,
    betalambda = betalambda,
    Z = Z
  )
  
  # 結果の出力 ---------------------------------
  return(res)
}

# 限界効果とデルタ法による標準誤差計算する関数
heckitmfx_delta <- function(model, newdata = NULL, type = c("unconditional", "conditional", "selection"), ...){
  
  type <- rlang::arg_match(type)
  # 限界効果の推定
  estimate <- heckitmfx(model = model, newdata = newdata, type = type) |> 
    select(term, all_of(type))
  
  # ヤコビ行列の計算
  grad <- numDeriv::jacobian(
    func = \(params) heckitmfx(.params = params, model = model, type = type, newdata = newdata) |> 
      pull(type), 
    x = coef(model) |> as_vector()
  )   
  if(checkmate::testClass(model, "heckitrob")){
    # モデルが heckitrob クラスなら第1段階と第2段階の共分散行列を結合し
    # モデル全体の共分散行列を作成します。
    vcv1 <- vcov(model$stage1)
    vcv2 <- vcov(model$stage2)
    
    O <- matrix(0, ncol = ncol(vcv2), nrow = nrow(vcv1))
    vcv_param <- rbind(cbind(vcv1, O), cbind(t(O), vcv2))
  }else{
    vcv_param <- vcov(model)
    # `sampleSelection::heckit` なら首座対角行列以外の部分が NA になるため0で置換
    vcv_param[is.na(vcv_param)] <- 0
  }
  
  # 疑問が残りますが、なぜか上手く(?)計算できています。
  # デルタ法による標準誤差の計算
  se <- sqrt(diag(grad %*% vcv_param %*% t(grad)))
  
  # 結果の出力
  res <-
    as_tibble(estimate) |>
    select(term, estimate = all_of(type)) |>
    mutate(
      std.error = se,
      statistic = estimate / std.error,                  # Z統計量
      p.value = 2*pnorm(abs(statistic), lower.tail = F) # 両側P値
    ) |>
    mutate(type = type, .after = term)
  return(res)
}