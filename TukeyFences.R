# compute_stats_from_freq -- median is defined as the 50th percentile of the expanded data
# (i.e., counts are taken into account). Input: data.frame with columns 'value' and 'count'.
#
# Args:
# - df: data.frame with numeric 'value' and numeric/integer 'count'
# - probs: vector of probabilities for weighted quantiles (default includes .25, .5, .75)
# - interp: if TRUE, interpolate inside a bin; if FALSE, use "lower" style (first value with cumulative >= target)
# - na.rm: remove NA rows if TRUE
#
# Returns:
# - total_count, mean, modes, median (weighted 50th percentile), weighted_quantiles, IQR (weighted),
#   unique_value_count, aggregated (value/count table)

compute_stats_from_freq <- function(df,
                                    probs = c(0.25, 0.5, 0.75),
                                    interp = TRUE,
                                    na.rm = TRUE) {
  if (!all(c("value", "count") %in% names(df))) {
    stop("Data frame must contain columns named 'value' and 'count'.")
  }
  
  v <- df$value
  cts <- df$count
  
  if (!is.numeric(v) || !is.numeric(cts)) stop("'value' and 'count' must be numeric.")
  if (any(is.na(v)) || any(is.na(cts))) {
    if (na.rm) {
      ok <- !(is.na(v) | is.na(cts))
      v <- v[ok]; cts <- cts[ok]
    } else stop("NA present in 'value' or 'count'. Set na.rm=TRUE to remove.")
  }
  if (any(cts < 0)) stop("'count' must be non-negative.")
  
  # Aggregate identical values (if needed) and sort by value
  agg <- aggregate(cts, by = list(value = v), FUN = sum)
  names(agg)[2] <- "count"
  ord <- order(agg$value)
  values <- agg$value[ord]
  counts <- agg$count[ord]
  
  total_count <- sum(as.numeric(counts))
  if (total_count == 0) stop("Total count is zero; cannot compute statistics.")
  
  # Weighted mean
  mean_val <- sum(as.numeric(values) * as.numeric(counts)) / total_count
  
  # Mode(s) by count
  max_count <- max(counts)
  modes <- values[counts == max_count]
  
  # Cumulative counts for weighted quantiles
  cum_counts <- cumsum(as.numeric(counts))
  
  # Weighted quantile (single probability p in [0,1])
  weighted_quantile_single <- function(p) {
    if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) stop("p must be in [0,1].")
    # target position in expanded data
    target <- p * total_count
    # special cases
    if (target <= 0) return(values[1])
    if (target >= total_count) return(tail(values, 1))
    idx <- which(cum_counts >= target)[1]
    if (is.na(idx)) return(tail(values, 1))
    if (!interp) {
      return(values[idx])
    } else {
      c_prev <- if (idx > 1) cum_counts[idx - 1] else 0
      # if this bin has zero count or target equals previous cumulative -> return current value
      if (counts[idx] == 0 || target == c_prev) return(values[idx])
      prop <- (target - c_prev) / counts[idx]  # fraction inside this bin
      v_prev <- if (idx > 1) values[idx - 1] else values[1]
      v_curr <- values[idx]
      return(v_prev + prop * (v_curr - v_prev))
    }
  }
  
  # Ensure 0.5 is included so median is the weighted 50th percentile
  probs_unique <- unique(c(probs, 0.5))
  qs <- vapply(probs_unique, weighted_quantile_single, numeric(1), USE.NAMES = FALSE)
  names(qs) <- as.character(probs_unique)
  
  weighted_median <- as.numeric(qs[as.character(0.5)])
  
  # Weighted IQR using 0.25 and 0.75
  q1 <- if ("0.25" %in% names(qs)) qs["0.25"] else weighted_quantile_single(0.25)
  q3 <- if ("0.75" %in% names(qs)) qs["0.75"] else weighted_quantile_single(0.75)
  iqr_weighted <- as.numeric(q3 - q1)
  
  list(
    total_count = total_count,
    mean = mean_val,
    modes = modes,
    median = weighted_median,             # median = 50th percentile considering counts
    weighted_quantiles = setNames(as.numeric(qs), names(qs)),
    IQR = iqr_weighted,
    unique_value_count = length(values),
    aggregated = data.frame(value = values, count = counts)
  )
}

# Example usage:
df <- data.frame(value = c(1,2,3,5,10), count = c(1e9, 2e9, 1e9, 5e9, 32e9))
res <- compute_stats_from_freq(df)
res$median   # weighted median (50th percentile of expanded data)
res$mean
res$IQR


library(dplyr)
df <- df %>%
  mutate(count_formatted = format(count, big.mark = ",", scientific = FALSE))
