# Parsing of CCI identifiers.
#
# Current format:  "Sending_subclass:Ligand->Receiving_subclass:Receptor"
# Legacy format:   "Sending_subclass:Ligand^Receiving_subclass:Receptor",
#                  optionally prefixed with "L_" (sending) or "R_" (receiving).

DIRECTIONS <- c("sending", "receiving")

#' Split CCI IDs into sending/receiving subclass, ligand and receptor.
parse_cci_ids <- function(ids) {
  ids <- as.character(ids)
  clean <- sub("^[LR]_", "", ids)
  parts <- strsplit(clean, "->|\\^")
  sides_ok <- lengths(parts) == 2L
  if (any(!sides_ok)) {
    stop("Unrecognised CCI IDs (expected 'Sub:Ligand->Sub:Receptor'): ",
         paste(head(ids[!sides_ok]), collapse = ", "))
  }
  send <- vapply(parts, `[`, character(1), 1L)
  recv <- vapply(parts, `[`, character(1), 2L)
  colon_ok <- grepl(":", send, fixed = TRUE) & grepl(":", recv, fixed = TRUE)
  if (any(!colon_ok)) {
    stop("CCI IDs without 'Subclass:Gene' on both sides: ",
         paste(head(ids[!colon_ok]), collapse = ", "))
  }
  data.table(
    cci = ids,
    sending_subclass = sub(":.*$", "", send),
    ligand = sub("^[^:]*:", "", send),
    receiving_subclass = sub(":.*$", "", recv),
    receptor = sub("^[^:]*:", "", recv)
  )
}

#' Parse CCI IDs and add the "anchor" for a QTL direction: the side of the
#' interaction whose gene coordinates the QTL was mapped at (ligand for
#' sending, receptor for receiving), plus class/superclass annotations.
annotate_cci <- function(ids, direction, ct) {
  direction <- match.arg(direction, DIRECTIONS)
  x <- parse_cci_ids(ids)
  x[, direction := direction]
  if (direction == "sending") {
    x[, `:=`(subclass = sending_subclass, gene = ligand,
             partner_subclass = receiving_subclass, partner_gene = receptor)]
  } else {
    x[, `:=`(subclass = receiving_subclass, gene = receptor,
             partner_subclass = sending_subclass, partner_gene = ligand)]
  }
  x[, `:=`(
    class = map_subclass(subclass, ct, "class"),
    superclass = map_subclass(subclass, ct, "superclass"),
    sending_class = map_subclass(sending_subclass, ct, "class"),
    receiving_class = map_subclass(receiving_subclass, ct, "class"),
    lr = paste0(ligand, "->", receptor)
  )]
  x[, cell_pair := paste0(sending_class, "->", receiving_class)]
  x[]
}
