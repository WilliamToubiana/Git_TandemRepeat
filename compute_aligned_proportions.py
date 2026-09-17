from pathlib import Path
import pybedtools
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# Set directories
WGADIR = Path("/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/WGA")
WINDOWDIR = Path("/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/windows")
ANNOTATIONDIR = Path("/work/FAC/FBM/DEE/tschwand/default/jsouto/will_TR_paper/annotation")

# Load full BEDs once (they're large and reused)
tr_coords = pybedtools.BedTool(ANNOTATIONDIR / "Tps_Tr.per250000w.sorted.bed")
aligned_coords = pybedtools.BedTool(WGADIR / "TPS_TO_TCM.TR.per250000w.sorted.n_aligned.bed")

# Convert to in-memory lists (safe for multiprocessing)
tr_coords_list = list(tr_coords)
aligned_coords_list = list(aligned_coords)

def process_window(window_line):
    win_chr = window_line.chrom
    win_start = window_line.start
    win_end = window_line.end

    # Total Tr bases overlapping this window
    total_tr = sum([
        min(f.end, win_end) - max(f.start, win_start)
        for f in tr_coords_list
        if f.chrom == win_chr and f.end > win_start and f.start < win_end
    ])

    # Aligned bases — now correctly calculated from 3rd - 2nd column
    aligned_bases = sum([
        min(f.end, win_end) - max(f.start, win_start)
        for f in aligned_coords_list
        if f.chrom == win_chr and f.end > win_start and f.start < win_end
    ])

    prop = aligned_bases / total_tr if total_tr > 0 else 0
    return [win_chr, win_start, win_end, aligned_bases, total_tr, round(prop, 5)]



# Load windows
windows = list(pybedtools.BedTool(WINDOWDIR / "Tps_w250000.sorted.bed"))

# Run with parallel processing using 4 workers
with ProcessPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(process_window, windows))

# Write results
df = pd.DataFrame(results, columns=["chrom", "start", "end", "aligned", "total", "prop"])
df.to_csv("proportion_TR_aligned_per_window.bed", sep="\t", index=False)
