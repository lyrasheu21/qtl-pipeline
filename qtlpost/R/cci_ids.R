# CCI identifiers, used by the CCI phenotype adapter and CCI-specific plots.
#
# Current format:  "Sending_subclass:Ligand->Receiving_subclass:Receptor"
# Legacy format:   "Sending_subclass:Ligand^Receiving_subclass:Receptor",
#                  optionally prefixed with "L_" (sending) or "R_" (receiving).

CCI_DIRECTIONS <- c("sending", "receiving")

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
    stop("CCI IDs without 'Subclass:Gene' on both sides: ", paste(head(ids[!colon_ok]), collapse = ", "))
  }
  data.table(
    cci = ids,
    sending_group = sub(":.*$", "", send),
    ligand = sub("^[^:]*:", "", send),
    receiving_group = sub(":.*$", "", recv),
    receptor = sub("^[^:]*:", "", recv)
  )
}

#' Phenotype annotation for CCIs mapped in one direction. The anchor (`gene`,
#' `group`) is the side whose gene coordinates the QTL is mapped at: the ligand
#' for sending, the receptor for receiving; `partner_*` is the other side.
cci_annotation <- function(ids, direction) {
  direction <- match.arg(direction, CCI_DIRECTIONS)
  x <- parse_cci_ids(ids)
  if (direction == "sending") {
    x[, `:=`(group = sending_group, gene = ligand, partner_group = receiving_group, partner_gene = receptor)]
  } else {
    x[, `:=`(group = receiving_group, gene = receptor, partner_group = sending_group, partner_gene = ligand)]
  }
  setnames(x, "cci", "phenotype_id")
  x[]
}
