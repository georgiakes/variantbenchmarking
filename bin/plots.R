#!/usr/bin/env Rscript

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

# Load necessary libraries
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(scales))
suppressWarnings(library(dplyr))

# Function to generate plots
generate_plots <- function(table, benchmark, type, filter, stats) {
    if (type != "None" && filter != "None") {
        table <- table[table$Type == type & table$Filter == filter, ]
        name1 <- paste(type, "_", filter, "_f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste(type, "_", filter, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste(type, "_", filter, "_pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    } else if (stats != "None") {
        table <- table[table$StatsType == stats, ]
        name1 <- paste(stats, "_f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste(stats, "_variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste(stats, "_pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    } else {
        name1 <- paste("f1_by_tool_", benchmark, "_mqc.png", sep = "")
        name2 <- paste("variants_by_tool_", benchmark, "_mqc.png", sep = "")
        name3 <- paste("pr_recall_by_tool_", benchmark, "_mqc.png", sep = "")
    }

    # Ensure Tool is a factor with consistent levels across all plots
    table$Tool <- factor(table$Tool)

    input_data_melted <- melt(table, id.vars = "Tool")

    tp_data <- input_data_melted[input_data_melted$variable %in% c( "TP_comp", "FP", "FN"), ]
    metric_data <- input_data_melted[input_data_melted$variable %in% c("F1"), ]
    metric_data$value <- as.numeric(as.character(metric_data$value))
    tp_data$value <- as.numeric(as.character(tp_data$value))

    # Specify the order of levels for the variable aesthetic
    tp_data$variable <- factor(tp_data$variable, levels = c("TP_comp", "FP", "FN"))
    metric_data$variable <- factor(metric_data$variable, levels = c("F1"))

    axis_text_size <- 12
    axis_title_size <- 14
    point_size <- 3
    line_size <- 1.2
    title_size <- 16
    facet_text_size <- 14
    legend_text_size <- 12
    legend_title_size <- 14

    # Calculate dynamic ranges for all plots
    # F1 plot range - always between 0 and 1
    f1_range <- range(metric_data$value, na.rm = TRUE)
    if (diff(f1_range) < 0.05) {
        f1_limits <- c(max(0, f1_range[1] - 0.02), min(1, f1_range[2] + 0.02))
    } else {
        f1_buffer <- max((f1_range[2] - f1_range[1]) * 0.1, 0.01)
        f1_limits <- c(max(0, f1_range[1] - f1_buffer), min(1, f1_range[2] + f1_buffer))
    }

    # Precision-Recall plot range - always between 0 and 1
    pr_range_x <- range(table$Recall, na.rm = TRUE)
    pr_range_y <- range(table$Precision, na.rm = TRUE)

    if (diff(pr_range_x) < 0.05) {
        pr_limits_x <- c(max(0, pr_range_x[1] - 0.02), min(1, pr_range_x[2] + 0.02))
    } else {
        pr_buffer_x <- max((pr_range_x[2] - pr_range_x[1]) * 0.1, 0.01)
        pr_limits_x <- c(max(0, pr_range_x[1] - pr_buffer_x), min(1, pr_range_x[2] + pr_buffer_x))
    }

    if (diff(pr_range_y) < 0.05) {
        pr_limits_y <- c(max(0, pr_range_y[1] - 0.02), min(1, pr_range_y[2] + 0.02))
    } else {
        pr_buffer_y <- max((pr_range_y[2] - pr_range_y[1]) * 0.1, 0.01)
        pr_limits_y <- c(max(0, pr_range_y[1] - pr_buffer_y), min(1, pr_range_y[2] + pr_buffer_y))
    }

    # Custom scale function for facets
    tp_plot <- ggplot(tp_data, aes(x = Tool, y = value, color = Tool)) +
    geom_line(aes(group = Tool), linewidth = line_size) +
    geom_point(size = point_size) +
    labs(x = "Tool", y = "Value", color = "Tool", title = "Variant Comparison Metrics") +
    facet_wrap(~variable, scales = "free_y") +
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size, margin = margin(t = 5)),
        axis.text.y = element_text(size = axis_text_size),
        axis.title.x = element_text(size = axis_title_size, margin = margin(t = 10)),
        axis.title.y = element_text(size = axis_title_size),
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        strip.text = element_text(size = facet_text_size, face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
        panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(
        labels = scales::label_number(big.mark = ",", accuracy = 1)
    )

    # Visualize f1 with clear y-axis values
    f1_plot <- ggplot(metric_data, aes(x = Tool, y = value, color = Tool)) +
        geom_point(size = point_size) +
        labs(x = "Tool", y = "F1 Score", title = "F1 Score") +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size, margin = margin(t = 5)),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size, margin = margin(t = 10)),
            axis.title.y = element_text(size = axis_title_size),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
            panel.grid.minor = element_blank()) +
        scale_y_continuous(
            labels = scales::label_number(accuracy = 0.001),
            limits = f1_limits,
            breaks = scales::breaks_extended(n = 8)
        ) +
        scale_color_discrete()

    # Visualize Precision vs Recall with clear axis values
    pr_plot <- ggplot(table) +
        geom_point(aes(x = Recall, y = Precision, color = Tool), size = point_size) +
        labs(x = "Recall", y = "Precision", title = "Precision vs Recall") +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size),
            axis.title.x = element_text(size = axis_title_size),
            axis.title.y = element_text(size = axis_title_size),
            plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_title_size),
            panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
            panel.grid.minor = element_blank()
        ) +
        scale_x_continuous(
            limits = pr_limits_x,
            labels = scales::label_number(accuracy = 0.001),
            breaks = scales::breaks_extended(n = 8)
        ) +
        scale_y_continuous(
            limits = pr_limits_y,
            labels = scales::label_number(accuracy = 0.001),
            breaks = scales::breaks_extended(n = 8)
        )

    # Save the plots
    tryCatch({
        if (!is.null(f1_plot)) {
            ggsave(name1, f1_plot, width = 6, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving metric plot: ", conditionMessage(e))
    })

    tryCatch({
        if (!is.null(tp_plot)) {
            ggsave(name2, tp_plot, width = 12, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving TP plot: ", conditionMessage(e))
    })

    tryCatch({
        if (!is.null(pr_plot)) {
            ggsave(name3, pr_plot, width = 6, height = 6, units = "in", dpi = 300, limitsize = TRUE)
        }
    }, error = function(e) {
        message("Error occurred while saving precision recall plot: ", conditionMessage(e))
    })
}

# Main script
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Please input the .summary file and benchmark", call. = FALSE)
}

table <- read.csv(args[1])
benchmark <- args[2]

if (benchmark == "happy") {
    generate_plots(table, benchmark, "SNP", "PASS", "None")
    generate_plots(table, benchmark, "SNP", "ALL", "None")
    generate_plots(table, benchmark, "INDEL", "PASS", "None")
    generate_plots(table, benchmark, "INDEL", "ALL", "None")
}else if (benchmark == "wittyer") {
    generate_plots(table, benchmark, "None", "None", "Base")
    generate_plots(table, benchmark, "None", "None", "Event")

}else {
    if (benchmark == "rtgtools") {
        table <- table[table$Threshold == "None", ]
    }
    generate_plots(table, benchmark, "None", "None", "None")
}
