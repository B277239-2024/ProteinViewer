# === evaluation/viz_minimal.R ===
# 读取 evaluation/ 内最新的 shiny.tictoc CSV
# 仅输出三张图：Top10 组件总耗时、事件耗时直方图、步骤内延迟（Per-step latency）

## ---------- 路径与输出目录 ----------
base_dir <- "evaluation"
fig_dir  <- file.path(base_dir, "figs")
tab_dir  <- file.path(base_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## 自动选择 evaluation/ 下最新的 csv（也可手动指定 csv_file <- "evaluation/xxx.csv"）
csv_candidates <- list.files(base_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_candidates) == 0) stop("evaluation/ 下找不到 CSV 文件。")
csv_file <- csv_candidates[which.max(file.info(csv_candidates)$mtime)]
message("使用 CSV：", csv_file)

## ---------- 依赖 ----------
need <- c("readr","dplyr","stringr","ggplot2","forcats","scales")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

## ---------- 读取与规范列 ----------
df <- readr::read_csv(csv_file, show_col_types = FALSE)

# 兼容列名（期望：measurement_id, type, start_time, duration (ms)）
nm <- names(df)
if (!"measurement_id" %in% nm && "name" %in% nm) {
  names(df)[match("name", nm)] <- "measurement_id"
}
if (!"duration (ms)" %in% nm && "duration" %in% nm) {
  names(df)[match("duration", names(df))] <- "duration (ms)"
}
stopifnot(all(c("measurement_id","type","start_time","duration (ms)") %in% names(df)))

df <- df |>
  rename(
    name        = measurement_id,
    duration_ms = `duration (ms)`
  ) |>
  mutate(
    end_time  = start_time + duration_ms,
    # 组件名：取 measurement_id 在 ":" 或 "-" 之前的前缀
    component = stringr::str_replace(name, "[:\\-].*$", "")
  )

## ---------- 映射到“你的 10 步” ----------
assign_step <- function(nm) {
  dplyr::case_when(
    # 1. Basic Info
    stringr::str_detect(nm, regex("^(info|domain_table|custom_domain_table_view)$", TRUE)) ~ 1L,
    # 2. gnomAD
    stringr::str_detect(nm, regex("^gnomad_(summary|table)$", TRUE)) ~ 2L,
    # 3. ConSurf（上传/图）
    stringr::str_detect(nm, regex("^(consurf1-|consurf_plot_ui|consurf1-txt_plot|consurf$)", TRUE)) ~ 3L,
    # 4. AlphaMissense
    stringr::str_detect(nm, regex("^(alphamissense_heatmap$|am1-|am_score_bar$)", TRUE)) ~ 4L,
    # 5. FELLS
    stringr::str_detect(nm, regex("^(fells_button$|download_fells_plot$|fells_plot$)", TRUE)) ~ 5L,
    # 6. 1D Plot
    stringr::str_detect(nm, regex("^(variant_1dplot$|plot_legend$|download_1dplot$)", TRUE)) ~ 6L,
    # 7. 3D 初始化
    stringr::str_detect(nm, regex("^(structure_view$|r3dmol:setSlab$)", TRUE)) ~ 7L,
    # 8. Sphere（注意：这里**不要**包含 addStyle）
    stringr::str_detect(nm, regex("^(r3dmol:addSphere$|sphere_legend$)", TRUE)) ~ 8L,
    # 9. ConSurf → 3D 着色（**把 addStyle / setStyle 都放到这里**）
    stringr::str_detect(nm, regex("^r3dmol:(addStyle|setStyle).*|^toggle_consurf_3d$", TRUE)) ~ 9L,
    # 10. Comparison
    stringr::str_detect(nm, regex("^comparison1-", TRUE)) ~ 10L,
    TRUE ~ NA_integer_
  )
}
df <- df |> mutate(step = assign_step(name))

## =========================================================
## 图 1：Top 10 components by total time  →  top10_components_total_time.png
## =========================================================
by_component <- df |>
  group_by(component) |>
  summarise(total_ms = sum(duration_ms, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(total_ms))

pA <- by_component |>
  slice_head(n = 10) |>
  mutate(component = forcats::fct_reorder(component, total_ms)) |>
  ggplot(aes(x = component, y = total_ms/1000)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Total duration (s)", title = "Top 10 components by total time")
ggsave(file.path(fig_dir, "top10_components_total_time.png"),
       pA, width = 8, height = 5, dpi = 200)

# 可选：导出组件表（便于在文中引用数字）
readr::write_csv(by_component, file.path(tab_dir, "total_time_by_component.csv"))

## =========================================================
## 图 2：Distribution of event durations  →  duration_histogram.png
## =========================================================
p50 <- stats::median(df$duration_ms, na.rm=TRUE)
p90 <- stats::quantile(df$duration_ms, 0.9, na.rm=TRUE)

pB <- ggplot(df, aes(x = duration_ms)) +
  geom_histogram(bins = 80) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si(" "))) +
  geom_vline(xintercept = p50, linetype = 2) +
  geom_vline(xintercept = p90, linetype = 3) +
  annotate("text", x = p50, y = Inf, vjust = 1.5, size = 3,
           label = sprintf("Median = %.0f ms", p50)) +
  annotate("text", x = p90, y = Inf, vjust = 3.0, size = 3,
           label = sprintf("P90 = %.0f ms", p90)) +
  labs(x = "Event duration (ms, log10)", y = "Count",
       title = "Distribution of event durations")
ggsave(file.path(fig_dir, "duration_histogram.png"),
       pB, width = 8, height = 5, dpi = 200)

## =========================================================
## 图 3：Per-step latency（步骤内延迟） →  per_step_latency.png
##   逻辑：每步 t0 = 该步内最早 start_time；
##        t1 = 该步“关键里程碑事件”的最早 end_time（找不到则用该步内最早 end_time）；
##        latency = (t1 - t0) / 1000 秒
## =========================================================
COMPLETION_MODE <- "last"  # 也可 "p95" 或 "first"

step_labels_short <- c(
  `1`="Basic", `2`="gNOMAD", `3`="ConSurf→1D",
  `4`="AM→1D", `5`="FELLS→1D",
  `6`="1D", `7`="3D Init",
  `8`="Sphere", `9`="ConSurf3D",
  `10`="Comp"
)

# 更严格、不会“越级提前”的锚点模式
anchor_patterns <- list(
  `1`  = c("^info$", "^domain_table$", "^custom_domain_table_view$"),
  `2`  = c("^gnomad_summary$", "^gnomad_table$"),
  `3`  = c("^consurf1-txt_plot$", "^consurf_plot_ui$", "^consurf$"),
  `4`  = c("^alphamissense_heatmap$", "^am1-", "^am_score_bar$"),
  `5`  = c("^fells_button$", "^download_fells_plot$", "^fells_plot$"),
  `6`  = c("^variant_1dplot$"),
  `7`  = c("^structure_view$", "^r3dmol:setSlab$"),
  `8`  = c("^r3dmol:addSphere$", "^sphere_legend", "sphere"),  # 放宽：包含 sphere 字样
  # 9 不再用 setStyle 做锚点（避免提前到 3D Init），只用 toggle_consurf_3d 或 addStyle
  `9`  = c("^toggle_consurf_3d$", "^r3dmol:addStyle"),
  `10` = c("^comparison1-")
)

# 顺序寻找锚点：第 k 步锚点必须 > 第 k-1 步锚点
find_anchor_after <- function(patterns, after_t) {
  hits <- do.call(rbind, lapply(patterns, function(pat){
    rows <- df[stringr::str_detect(df$name, pat) & df$start_time > after_t, c("start_time","end_time")]
    if (nrow(rows)) rows else NULL
  }))
  if (is.null(hits) || nrow(hits) == 0) return(NA_real_)
  min(hits$start_time, na.rm = TRUE)
}

steps_seq <- as.integer(names(step_labels_short))
anchors <- data.frame(step = steps_seq,
                      label = unname(step_labels_short),
                      start_t = NA_real_)

prev_t <- -Inf
for (i in seq_along(steps_seq)) {
  s <- steps_seq[i]
  t <- find_anchor_after(anchor_patterns[[as.character(s)]], prev_t)
  anchors$start_t[i] <- t
  if (!is.na(t)) prev_t <- t
}
# 去掉未找到锚点的步骤
anchors <- anchors[!is.na(anchors$start_t), ]
anchors <- anchors[order(anchors$start_t), , drop = FALSE]
anchors$end_t <- c(anchors$start_t[-1], max(df$end_time, na.rm = TRUE))

calc_completion <- function(sub_df, mode="last"){
  if (nrow(sub_df) == 0) return(NA_real_)
  if (mode == "last") {
    return(max(sub_df$end_time, na.rm = TRUE))
  } else if (mode == "p95") {
    ord <- sub_df[order(sub_df$end_time), ]
    cs  <- cumsum(ord$duration_ms); thr <- 0.95 * sum(ord$duration_ms)
    return(ord$end_time[which(cs >= thr)[1]])
  } else if (mode == "first") {
    return(min(sub_df$end_time, na.rm = TRUE))
  } else stop("Unknown COMPLETION_MODE")
}

latency_df <- do.call(rbind, lapply(seq_len(nrow(anchors)), function(i){
  t0 <- anchors$start_t[i]; t1 <- anchors$end_t[i]
  win_df <- df[df$start_time >= t0 & df$start_time < t1, , drop = FALSE]
  end_pt <- calc_completion(win_df, COMPLETION_MODE)
  data.frame(step = anchors$step[i],
             label = anchors$label[i],
             latency_s = (end_pt - t0)/1000)
}))

# 导出表 + 作图
readr::write_csv(latency_df, file.path(tab_dir, "per_step_latency.csv"))

pC2 <- ggplot2::ggplot(latency_df,
                       ggplot2::aes(x = factor(label, levels = step_labels_short), y = latency_s)) +
  ggplot2::geom_col() +
  ggplot2::labs(x = NULL, y = "Latency after step start (s)",
                title = "Per-step completion time") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
ggplot2::ggsave(file.path(fig_dir, "per_step_latency.png"),
                pC2, width = 10, height = 5, dpi = 200)

message("完成。图像已输出到：", normalizePath(fig_dir))
message("相关数据表输出到：", normalizePath(tab_dir))
