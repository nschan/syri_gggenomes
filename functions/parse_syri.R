#' Parse syri output
#' Parse syri output to plot with gggenomes / geom_polygon
#' 
#' @param files 
#' a list of files. These files are expected to: end with `.syri.out` and follow a naming scheme like A_on_B.syri.out
#' @param order 
#' a dataframe with a column bin_id , containing the order of genomes
#' @param resize_polygons 
#' (logical) should polygons of short links be resized?
#' @param resize_polygons_size 
#' if polygons are resized, to what fraction of the total length? Default `0.003`
#' @param min_polygon_feat_size 
#' minimum length of links to be resized
#' @param no_polygons
#' do not compute polygons (default: false, will compute polygons)
#'
#' @return 
#' a list of dataframes: 
#' `$seqs`: contains sequenece information, according to gggenomes
#' `$links`: contains links between sequences, according to gggenomes
#' `$polygons`: contains polygons that can be plotted via `geom_polygon()`, might be prettier than `links`
#' 
#' @export
parse_syri <- \(
  files,
  order,
  resize_polygons = T,
  resize_polygons_size = 0.003,
  min_polygon_feat_size = 5000,
  no_polygons = FALSE
) {
  seqs <- lapply(files, \(file) {
    seq_names <- file %>%
      stringr::str_remove('.+?(?=[A-Za-z0-9_-]*.syri.out)')  %>%
      stringr::str_remove_all(".syri.out") %>%
      stringr::str_split("_on_", simplify = T)
    
    ref <- seq_names[[2]]
    qry <- seq_names[[1]]
    #message("Got ref and query names")
    syri_file <- vroom::vroom(
      file,
      col_names = c(
        "seq_id",
        "ref_start",
        "ref_end",
        "ref_seq",
        "query_seq",
        "seq_id2",
        "query_start",
        "query_end",
        "unique_ID",
        "parent_ID",
        "annotation",
        "copy"
      ),
      show_col_types = FALSE
    ) %>%
      dplyr::filter(annotation %in% c("SYN", "INV", "TRANS", "DUP")) %>%
      dtplyr::lazy_dt() %>%
      dplyr::mutate(bin_id = ref, bin_id2 = qry) %>%
      dplyr::mutate(across(ends_with(c("start", "end")), as.integer)) %>%
      dplyr::mutate(
        length = max(ref_start, ref_end, na.rm = T),
        length2 = max(query_start, query_end, na.rm = T),
      )  %>%
      as.data.frame()
    #message("Formatted input")
    syri_head <- syri_file %>%
      slice_head(n = 1)
    seqtab <- tibble(
      seq_id = c(syri_head$seq_id, syri_head$seqid2),
      bin_id = c(ref, qry)
    ) %>%
      dplyr::left_join(bind_rows(
        dplyr::select(syri_head, bin_id , seq_id, length),
        dplyr::select(
          syri_head,
          bin_id = bin_id2,
          seq_id = seq_id2,
          length = length2
        )
      ),
      by = c("bin_id", "seq_id"))
    return(seqtab)
  }) %>%
    bind_rows()
  message("Created seqtab")
  
  links <- lapply(files, \(file) {
    seq_names <- file %>%
      stringr::str_remove('.+?(?=[A-Za-z0-9_-]*.syri.out)') %>%
      stringr::str_remove_all(".syri.out") %>%
      stringr::str_split("_on_", simplify = T)
    
    ref <- seq_names[[2]]
    qry <- seq_names[[1]]
    
    syri_in <- suppressWarnings(vroom::vroom(
      file,
      col_names = c(
        "seq_id",
        "ref_start",
        "ref_end",
        "ref_seq",
        "query_seq",
        "seq_id2",
        "query_start",
        "query_end",
        "unique_ID",
        "parent_ID",
        "annotation",
        "copy"
      ),
      show_col_types = FALSE
    )) %>%
      dtplyr::lazy_dt() %>%
      dplyr::filter(annotation %in% c("SYN", "INV", "TRANS", "DUP")) %>%
      dplyr::mutate(bin_id = ref, bin_id2 = qry) %>%
      dplyr::mutate(across(ends_with(c("start", "end")), as.integer)) %>%
      dplyr::mutate(
        length = max(ref_start, ref_end, na.rm = T),
        length2 = max(query_start, query_end, na.rm = T),
      ) %>%
      dplyr::select(
        bin_id,
        seq_id,
        length,
        start = ref_start,
        end = ref_end,
        bin_id2,
        seq_id2,
        length2,
        start2 = query_start,
        end2 = query_end,
        unique_ID,
        parent_ID,
        type = annotation
      ) %>%
      as.data.frame()
  }) %>%
    bind_rows()
  message("Created links")
  results <- list()
  results$seqs <- left_join(
    order,
    seqs %>%
      dplyr::group_by(seq_id, bin_id) %>%
      dplyr::summarize(length = max(length), .groups = "drop"),
    by = c("bin_id")
  )
  results$links <- links
  if(!no_polygons) {
    message("Calculating polygons")
    results$polygons <- compute_polygons_syri(
      plotdat = results,
      as.character(results$seqs$bin_id),
      resize_polygons = resize_polygons,
      resize_size = resize_polygons_size,
      min_feat_size = min_polygon_feat_size
    ) %>%
      mutate(type = fct_relevel(type, "SYN", "DUP", "TRANS", "INV"))
    results$links <- results$links %>%
      mutate(type = fct_relevel(type, "SYN", "DUP", "TRANS", "INV"))
  }
  return(results)
}

# Helper for compute polygons
get_genomes <- function(plotdat, genome_order) {
  genomes_numbers <-
    c(1:nrow(plotdat$seqs)) %>%
    as.list()
  names(genomes_numbers) <- rev(genome_order)
  genomes_numbers
}

#' Compute polygons
#' Compute polygons to draw links when plotting with gggenomes
#' Wraps GENESPACE::calc_curvePolygon.
#' @param plotdat list of dfs, named seqs (containing sequences) and links (will be converted to polygons)
#' @param genome_order a vector containing genomes in the order they should appear in the plot
#' @param resize_polygons bool, should small features be reiszed (to 0.3% of the chromosome length)
#' @param resize_size number, resize to which size? Default: 0.3% of the chromosome length)
#' @param min_feat_size mininum size (in bp) of features to be kept / resized, default: 5000
#'
#' @return a dataframe containing polygons
#'
#' @export
#'
#' @examples
compute_polygons_syri <- function(plotdat,
                                  genome_order,
                                  resize_polygons = T,
                                  resize_size = 0.003 ,
                                  min_feat_size = 5000) {
  # polygon resizing should make them at least  0.3% of the chromosome length
  seq_lengths <- plotdat$seqs
  # Get the order of genomes
  genomes <- get_genomes(plotdat, genome_order)
  # Compute all polygons
  
  parallel::mclapply(1:nrow(plotdat$links), \(row_num) {
  #lapply(1:nrow(plotdat$links), \(row_num) {
    tmpdat = plotdat$links[row_num, ]
    if (tmpdat$end - tmpdat$start > min_feat_size &
        tmpdat$end2 - tmpdat$start2 > min_feat_size) {
      mid_x = tmpdat %$% mean(c(start, start2, end, end2))
      mid_y = tmpdat %$% mean(c(genomes[[paste(bin_id)]], genomes[[paste(bin_id2)]]))
      min_len1 = seq_lengths %>% filter(bin_id == tmpdat$bin_id) %$% length * resize_size
      min_len2 = seq_lengths %>% filter(bin_id == tmpdat$bin_id2) %$% length * resize_size
      if (resize_polygons) {
        # Add 0.5% on either side.
        if (tmpdat$end - tmpdat$start < min_len1) {
          mid_1 = mean(c(tmpdat$start, tmpdat$end))
          tmpdat$start = mid_1 - min_len1 / 2
          tmpdat$end = mid_1 + min_len1 / 2
        }
        if (tmpdat$end2 - tmpdat$start2 < min_len2) {
          mid_2 = mean(c(tmpdat$start2, tmpdat$end2))
          tmpdat$start2 = mid_2 - min_len2 / 2
          tmpdat$end2 = mid_2 + min_len2 / 2
        }
      }
      # For inversions
      if (tmpdat$type == "INV") {
        polygons <- bind_rows(
          tmpdat %$%
            calc_curve_poly(
              start1 = start,
              end1 = end,
              start2 = mid_x - 1,
              end2 = mid_x + 1,
              y1 = genomes[[paste(bin_id)]],
              y2 = mid_y - 0.01,
              npts = 500
            ) %>%
            as.data.frame() %>%
            dplyr::mutate(
              link = paste(
                tmpdat$bin_id,
                tmpdat$seq_id,
                tmpdat$bin_id2,
                tmpdat$seq_id2,
                row_num
              ),
              link_grp = paste(
                tmpdat$bin_id,
                tmpdat$seq_id,
                tmpdat$bin_id2,
                tmpdat$seq_id2,
                row_num,
                "lower"
              ),
              type = "INV"
            ),
          tmpdat %$%
            calc_curve_poly(
              start1 = mid_x - 1,
              end1 = mid_x + 1,
              start2 = start2,
              end2 = end2,
              y1 = mid_y + 0.01,
              y2 = genomes[[paste(bin_id2)]],
              npts = 500
            ) %>%
            as.data.frame() %>%
            dplyr::mutate(
              link = paste(
                tmpdat$bin_id,
                tmpdat$seq_id,
                tmpdat$bin_id2,
                tmpdat$seq_id2,
                row_num
              ),
              link_grp = paste(
                tmpdat$bin_id,
                tmpdat$seq_id,
                tmpdat$bin_id2,
                tmpdat$seq_id2,
                row_num,
                "upper"
              ),
              type = "INV"
            )
        )
      } else {
        polygons <- tmpdat %$%
          calc_curve_poly(
            start1 = start,
            end1 = end,
            start2 = start2,
            end2 = end2,
            y1 = genomes[[paste(bin_id)]],
            y2 = genomes[[paste(bin_id2)]],
            npts = 1000
          ) %>%
          as_tibble() %>%
          dplyr::mutate(
            link = paste(
              tmpdat$bin_id,
              tmpdat$seq_id,
              tmpdat$bin_id2,
              tmpdat$seq_id2,
              row_num
            ),
            link_grp = paste(
              tmpdat$bin_id,
              tmpdat$seq_id,
              tmpdat$bin_id2,
              tmpdat$seq_id2,
              row_num
            ),
            type = tmpdat$type
          )
      }
    } else {
      polygons <- NA
    }
    return(polygons)
  }) %>%
    {
      .[!is.na(.)]
    } %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(link) %>%
    dplyr::mutate(direct = case_when(max(y) - min(y)  == 1 ~ TRUE, TRUE ~ FALSE)) %>%
    dplyr::ungroup()
}

calc_curve_poly <- function(start1,
                            end1 = NULL,
                            start2,
                            end2 = NULL,
                            y1,
                            y2,
                            npts = 250,
                            keepat = round(npts / 20)) {
  # This is GENESPACE::calc_curvePolygon, to not depend on GENESPACE package for this single function.
  cos_points <- function(npts, keepat) {
    # initial number of points
    # grid to keep always
    grid <- seq(from = 0,
                to = pi,
                length.out = npts) # grid
    x <- (1 - cos(grid)) / max((1 - cos(grid))) # scaled cosine
    y <- grid / max(grid) # scaled grid
    # calculate slope for each point
    x1 <- x[-1]
    y1 <- y[-1]
    x2 <- x[-length(x)]
    y2 <- y[-length(y)]
    s <-  (y1 - y2) / (x1 - x2)
    # choose points that capture changes in slope
    ds <- cumsum(abs(diff(s))) * 5
    wh <- c(1, which(!duplicated(round(ds))), length(x))
    wh2 <- c(wh, seq(
      from = 0,
      to = length(x),
      by = round(keepat)
    ))
    wh <- c(wh, wh2)[!duplicated(c(wh, wh2))]
    wh <- wh[order(wh)]
    return(cbind(x[wh], y[wh]))
  }
  
  scaledCurve <- cos_points(npts = npts, keepat = keepat)
  # print(scaledCurve)
  if (!is.null(end1) | !is.null(end2)) {
    sc1 <- scaledCurve[, 1]
    sc2 <- scaledCurve[, 2]
    
    tp <- rbind(
      start1 = data.table::data.table(x = start1,
                                      y = y1),
      poly1 = data.table::data.table(
        x = scale_betwn(x = sc1, min = start1, max = start2),
        y = scale_betwn(x = sc2, min = y1, max = y2)
      ),
      start2 = data.table::data.table(x = start2, y = y2),
      end2 = data.table::data.table(x = end2, y = y2),
      poly2 = data.table::data.table(
        x = scale_betwn(x = sc1, min = end2, max = end1),
        y = scale_betwn(x = sc2, min = y2, max = y1)
      ),
      end1 = data.table::data.table(x = end1, y = y1)
    )
  } else{
    tp <- data.table::data.table(
      x = scale_betwn(x = scaledCurve[, 1], min = start1, max = start2),
      y = scale_betwn(x = scaledCurve[, 2], min = y1, max = y2)
    )
  }
  return(tp)
}

scale_betwn <- function(x, min, max, scale1toMean = TRUE) {
  # This is GENESPACE::scale_between, to not depend on GENESPACE package for this single function.
  if (length(unique(x)) > 1) {
    return((x - min(x)) / (max(x) - min(x)) * (max - min) + min)
  } else{
    if (scale1toMean) {
      return(mean(c(min, max)))
    } else{
      return(max)
    }
  }
}

# Plot colors
syri_plot_fills <- scale_fill_manual(
  labels = c(
    "SYN" = "Syntenic",
    "INV" = "Inversion",
    "TRANS" = "Translocation",
    "DUP" = "Duplication"
  ),
  values = c(
    "SYN" = "#bcbcbc",
    "INV" = "#7c02d3",
    "TRANS" = "#3dccf7",
    "DUP" = "#f9c131"
  ),
  name = "Relationship"
) 