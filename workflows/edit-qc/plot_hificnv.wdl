version 1.0
import "../wdl-common/wdl/structs.wdl"
task plot_CNV {
    input {
        File depth_bw
        File maf_bw
        Int bins_genome = 2000
        Int bins_chrom = 500 
        RuntimeAttributes runtime_attributes
    }

    Int threads   = 16
    Int mem_gb    = threads * 5
    Int disk_size = ceil(size(depth_bw, "GB") + size(maf_bw, "GB") + 100)


    command <<<
        set -euo pipefail
        mkdir -p tmp
        pip install --target=tmp matplotlib
        export PYTHONPATH="tmp"
        cat << EOF > plot_coverage.py 
        import re
        import pyBigWig
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mtick
        from collections import defaultdict

        # ─── Helper: natural contig ordering ───────────────────────────────────────────
        def order_contigs(chrom_lengths):
            grouped = defaultdict(list)
            for chrom in chrom_lengths:
                m = re.match(r'^(chr[^_]+)', chrom)
                base = m.group(1) if m else chrom
                grouped[base].append(chrom)

            def base_key(b):
                nb = b[3:] if b.startswith("chr") else b
                if   nb.isdigit():             return (0, int(nb))
                elif nb.upper() == "X":        return (1, 23)
                elif nb.upper() == "Y":        return (1, 24)
                elif nb.upper() in ("M","MT"): return (1, 25)
                else:                           return (2, b)

            ordered = []
            for base in sorted(grouped, key=base_key):
                parts = sorted(grouped[base], key=lambda x: (0 if x == base else 1, x))
                ordered.extend(parts)
            return ordered

        # ─── Helper: compute per‐contig binned coverage ────────────────────────────────
        def compute_profile(bw, contigs, bins=2000):
            vals = []
            bounds = [0]
            offset = 0
            sizes = bw.chroms()

            for c in contigs:
                L = sizes[c]
                s = bw.stats(c, 0, L, type="mean", nBins=bins) or [0]*bins
                vals.extend([0 if v is None else v for v in s])
                offset += L
                bounds.append(offset)

            coords = np.linspace(0, bounds[-1], len(vals), endpoint=False)
            return coords, vals, bounds

        # ─── load & prep ─────────────────────────────────────────────────────────────
        bw      = pyBigWig.open("~{depth_bw}")
        ordered = order_contigs(bw.chroms())
        primary = [c for c in ordered if re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2]|X|Y|M)", c)]
        coords, vals, bounds = compute_profile(bw, primary, bins=2000)
        mean_cov = np.mean(vals)
        # ─── plot ──────────────────────────────────────────────────────────────────────
        plt.rcParams.update({
            "figure.figsize": (20, 4),
            "axes.titlesize": 18,
            "axes.labelsize": 14,
            "ytick.labelsize": 10,
            "legend.fontsize": 12,
        })

        fig, ax = plt.subplots()
        ax.plot(coords, vals, lw=0.6)

        # remove horizontal padding so plot starts exactly at chr1
        ax.margins(x=0)
        # alternate shading by chromosome
        for i, (s, e) in enumerate(zip(bounds[:-1], bounds[1:])):
            if i % 2 == 0:
                ax.axvspan(s, e, color="gray", alpha=0.1)

        # fixed y‐limit + formatting
        ax.set_ylim(0, 60)
        ax.yaxis.set_major_locator(mtick.MultipleLocator(10))
        ax.yaxis.set_major_formatter(lambda x, pos: f"{int(x)}×")
        # mean line
        ax.axhline(mean_cov, color="red", ls="--")

        # remove default x‐ticks & labels
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

        # make room for custom ticks, labels, and legend
        plt.subplots_adjust(bottom=0.30, right=0.75)

        # custom chromosome‐boundary ticks
        ax.set_xticks(bounds)
        ax.tick_params(
            axis='x', 
            which='major', 
            bottom=True, 
            length=12,     # longer ticks
            width=2        # thicker ticks
        )

        # chromosome names centered in each interval
        label_y = -0.08  # 8% below axis in axes-fraction
        for start, end, name in zip(bounds[:-1], bounds[1:], primary):
            ax.text(
                x=(start + end) / 2,
                y=label_y,
                s=name,
                ha="center",
                va="top",
                rotation=45,
                fontweight='bold',
                transform=ax.get_xaxis_transform()
            )

        # axis labels
        ax.set_xlabel("Chromosome", labelpad=50)   # push xlabel below tick labels
        ax.set_ylabel("Average per-base coverage")

        # legend outside
        ax.set_title("Genome-wide profile (primary chromosomes only)")
        plt.tight_layout()
        fig.savefig("genome_profile.png", dpi=300)

        # per-chromosome (global-cords)
        import numpy as np
        import matplotlib.ticker as mtick

        for chrom in primary:
            # 1) get coords, values & boundaries
            c_coords, c_vals, c_bounds = compute_profile(bw, [chrom], bins=500)

            # 2) start figure
            fig, ax = plt.subplots(figsize=(12, 3))

            # 3) plot coverage
            ax.plot(c_coords, c_vals, lw=0.6)

            # 4) same y-range + formatting
            ax.set_ylim(0, 60)
            ax.yaxis.set_major_locator(mtick.MultipleLocator(10))
            ax.yaxis.set_major_formatter(lambda y, pos: f"{int(y)}×")

            # 5) full-chr x-range
            chrom_len = c_bounds[-1]
            ax.set_xlim(0, chrom_len)

            # 6) ~6 evenly spaced bp ticks
            ticks = np.linspace(0, chrom_len, 6)
            ax.set_xticks(ticks)
            ax.set_xticklabels(
                [f"{int(x):,}" for x in ticks],  # comma‐formatted
                rotation=45, ha="right", fontsize=10
            )
            ax.tick_params(axis='x', length=8, width=1.5, pad=8)

            # 7) push bottom up so labels are visible
            fig.subplots_adjust(bottom=0.30)

            # 8) labels & title
            ax.set_xlabel("Position (bp)", labelpad=12)
            ax.set_ylabel("Coverage")
            ax.set_title(chrom)
            fig.savefig(f"{chrom}.png", dpi=150)
        bw.close()

        ## MAF PLOT
        # ─── 1) Load & sample via stats() ─────────────────────────────────────────────
        bw = pyBigWig.open("~{maf_bw}")
        all_maf = []

        for chrom, length in bw.chroms().items():
            # sample every chromosome into 20k equal bins
            stats = bw.stats(chrom, 0, length, type="mean", nBins=20_000) or []
            # filter out None
            cleaned = [v for v in stats if v is not None]
            all_maf.append(cleaned)

        bw.close()
        maf_vals = np.concatenate(all_maf)

        # how many true zeros?
        n_zero  = np.sum(maf_vals == 0)
        n_total = len(maf_vals)
        print(f"{n_zero} / {n_total} bins ({100*n_zero/n_total:.1f}%) have MAF == 0")

        # ─── 2) Prepare the figure ────────────────────────────────────────────────────
        plt.style.use("ggplot")
        fig, ax = plt.subplots(figsize=(8, 5))

        # ─── 3) Plot on a log-y scale so the tail is visible ─────────────────────────
        ax.hist(
            maf_vals,
            bins=100,
            range=(0, 1),
            edgecolor="black",
            linewidth=1,
            alpha=0.8,
            log=True
        )

        # ─── 4) Axis formatting ───────────────────────────────────────────────────────
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 11))
        ax.set_xlabel("Minor allele frequency", fontsize=12)
        ax.set_ylabel("Count (log scale)", fontsize=12)
        ax.tick_params(axis="both", which="major", labelsize=10)

        # ─── 5) Title & annotation for zeros ─────────────────────────────────────────
        ax.set_title("Genome-wide MAF distribution", fontsize=14, pad=12)

        # annotate the zero-bin percentage
        ax.text(
            0.02, 0.95,
            f"{100*n_zero/n_total:.1f}% zeros",
            transform=ax.transAxes,
            va="top",
            fontsize=11,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.7)
        )

        # ─── 6) Grid and layout ──────────────────────────────────────────────────────
        ax.grid(axis="y", alpha=0.4)
        plt.tight_layout()
        fig.savefig("maf_distribution.png", dpi=300)
        bw.close()
        EOF

        python plot_coverage.py
    >>>

    output {
        File genome_profile = "genome_profile.png"
        File maf_distribution = "maf_distribution.png"
        Array[File] chrom_profiles = glob("chr*.png")

    }

    runtime {
        docker: "quay.io/biocontainers/pybigwig:0.3.24--py311hd8c7dd8_0"
        cpu: threads
        memory: mem_gb + " GB"
        disk: disk_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: runtime_attributes.preemptible_tries
        maxRetries: runtime_attributes.max_retries
        awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
        zones: runtime_attributes.zones

    }
}
