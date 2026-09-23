"""Pick the donors to drop for a donor-fraction sensitivity variant (seeded)."""

import random

donors = [line.strip() for line in open(snakemake.input[0]) if line.strip()]
n_drop = round(len(donors) * (1 - snakemake.params.keep_fraction))
rng = random.Random(f"{snakemake.params.seed}-{snakemake.wildcards.variant}")
dropped = sorted(rng.sample(donors, n_drop))
with open(snakemake.output[0], "w") as fh:
    fh.writelines(f"{d}\n" for d in dropped)
with open(snakemake.log[0], "w") as fh:
    fh.write(f"{len(donors)} donors, dropping {n_drop}\n")
