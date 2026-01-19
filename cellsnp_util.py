# Function to read the vcf file from cellSNP output and generate the coverage_df as well as the mutrate_df
def read_cellsnp_vcf_to_matrices(vcf_gz_file, chunk_size=100000, progress_interval=100000):
    """
    Optimized version: Read cellSNP.cells.vcf.gz file and generate coverage and mutation rate matrices.

    Optimizations:
    - Uses NumPy arrays instead of Python lists for faster operations
    - Pre-allocates arrays when possible
    - Vectorized calculations
    - Processes data in chunks to reduce memory overhead
    - Minimizes string operations and conversions

    Parameters:
    -----------
    vcf_gz_file : str
        Path to cellSNP.cells.vcf.gz file
    chunk_size : int, default=100000
        Number of positions to process in each chunk (helps with memory management)
    progress_interval : int, default=100000
        Print progress every N variants. Set to 0 to disable progress output.

    Returns:
    --------
    tuple: (coverage_df, mutrate_df)
        - coverage_df: DataFrame with total coverage
        - mutrate_df: DataFrame with mutation rates
    """
    import gzip
    import pandas as pd
    import numpy as np

    print(f"Reading VCF file: {vcf_gz_file}")

    col_names = []
    row_names = []
    coverage_chunks = []
    mutrate_chunks = []
    
    current_coverage = []
    current_mutrate = []
    n_cells = 0
    n_variants = 0

    with gzip.open(vcf_gz_file, 'rt') as f:
        for line in f:
            # Skip header lines starting with ##
            if line.startswith('##'):
                continue

            # Parse column names from the #CHROM line
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                format_idx = parts.index('FORMAT')
                col_names = parts[format_idx + 1:]
                n_cells = len(col_names)
                print(f"Found {n_cells} cells")
                continue

            # Parse data lines
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            row_name = f"{chrom}.{pos}"
            row_names.append(row_name)

            # Extract format fields once
            format_fields = parts[8].split(':')

            # Find indices for DP, AD, and ALL in the FORMAT string
            try:
                dp_idx = format_fields.index('DP')
                ad_idx = format_fields.index('AD')
                all_idx = format_fields.index('ALL')
            except ValueError as e:
                raise ValueError(f"Missing required FORMAT field at position {row_name}: {e}")

            # Pre-allocate arrays for this position
            coverage_row = np.zeros(n_cells, dtype=np.int32)
            mutrate_row = np.zeros(n_cells, dtype=np.float32)

            # Parse each cell's data (starting from column 9)
            for cell_idx, cell_data in enumerate(parts[9:]):
                # Handle missing data
                if cell_data.startswith('.'):
                    continue

                # Split cell data by ':'
                cell_fields = cell_data.split(':')

                try:
                    # Extract values
                    dp_str = cell_fields[dp_idx]
                    ad_str = cell_fields[ad_idx]
                    all_str = cell_fields[all_idx]
                    
                    # Fast path for missing values
                    if dp_str == '.' or all_str == '.':
                        continue
                    
                    dp_val = int(dp_str)
                    ad_val = int(ad_str) if ad_str != '.' else 0

                    # Parse ALL field efficiently: "A,C,G,T,N"
                    # Sum first 4 values (A,C,G,T) without creating list
                    all_parts = all_str.split(',')
                    coverage = sum(int(all_parts[i]) for i in range(min(4, len(all_parts))))

                    # Reference reads = DP - AD
                    ref_val = dp_val - ad_val

                    # Mutation = Coverage - Reference reads
                    mutation = coverage - ref_val

                    # Store values
                    coverage_row[cell_idx] = coverage
                    
                    # Mutation rate (vectorized division will be faster later)
                    if coverage > 0:
                        mutrate_row[cell_idx] = mutation / coverage

                except (ValueError, IndexError):
                    # Keep as 0
                    continue

            current_coverage.append(coverage_row)
            current_mutrate.append(mutrate_row)
            n_variants += 1

            # Progress tracking
            if progress_interval > 0 and n_variants % progress_interval == 0:
                print(f"  Processed {n_variants:,} variants...")

            # Process in chunks to manage memory
            if len(current_coverage) >= chunk_size:
                coverage_chunks.append(np.array(current_coverage))
                mutrate_chunks.append(np.array(current_mutrate))
                current_coverage = []
                current_mutrate = []

    # Add remaining data
    if current_coverage:
        coverage_chunks.append(np.array(current_coverage))
        mutrate_chunks.append(np.array(current_mutrate))

    print(f"Parsed {len(row_names)} positions")

    # Concatenate all chunks
    print("Creating DataFrames...")
    coverage_array = np.vstack(coverage_chunks) if len(coverage_chunks) > 1 else coverage_chunks[0]
    mutrate_array = np.vstack(mutrate_chunks) if len(mutrate_chunks) > 1 else mutrate_chunks[0]
    
    coverage_df = pd.DataFrame(coverage_array, index=row_names, columns=col_names)
    mutrate_df = pd.DataFrame(mutrate_array, index=row_names, columns=col_names)

    print("Processing complete!")
    print(f"Coverage matrix shape: {coverage_df.shape}")
    print(f"Mutation rate matrix shape: {mutrate_df.shape}")
    print(f"Coverage range: [{coverage_df.min().min()}, {coverage_df.max().max()}]")
    print(f"Mutation rate range: [{mutrate_df.min().min():.4f}, {mutrate_df.max().max():.4f}]")

    return coverage_df, mutrate_df

def read_cellsnp_vcf_to_matrices_pysam_sparse(vcf_gz_file, progress_interval=100000):
    """
    Memory-efficient reading of cellSNP.cells.vcf.gz using sparse matrices.

    This version stores only non-zero values, dramatically reducing memory
    usage for sparse single-cell data (typically >95% memory savings).

    Parameters:
    -----------
    vcf_gz_file : str
        Path to cellSNP.cells.vcf.gz file
    progress_interval : int, default=100000
        Print progress every N variants. Set to 0 to disable progress output.

    Returns:
    --------
    tuple: (coverage_sparse, mutrate_sparse, row_ids, col_names)
        - coverage_sparse: scipy.sparse.csr_matrix with coverage values
        - mutrate_sparse: scipy.sparse.csr_matrix with mutation rates
        - row_ids: list of position identifiers (chrom.pos)
        - col_names: list of cell/sample names
    """
    import pysam
    import numpy as np
    from scipy import sparse

    print(f"Reading VCF file (sparse mode): {vcf_gz_file}")

    vcf = pysam.VariantFile(vcf_gz_file)
    samples = list(vcf.header.samples)
    n_cells = len(samples)
    print(f"Found {n_cells} cells")

    row_ids = []
    cov_data, cov_rows, cov_cols = [], [], []
    mut_data, mut_rows, mut_cols = [], [], []

    for i, record in enumerate(vcf):
        row_ids.append(f"{record.chrom}.{record.pos}")

        for cell_idx, sample_name in enumerate(samples):
            sample = record.samples[sample_name]
            dp_val = sample.get('DP')
            ad_val = sample.get('AD')
            all_vals = sample.get('ALL')

            if dp_val is None or all_vals is None:
                continue

            try:
                if isinstance(all_vals, (list, tuple)) and len(all_vals) >= 4:
                    cov = sum(all_vals[:4])
                    if cov > 0:  # Only store non-zero values
                        ref = dp_val - (ad_val if ad_val else 0)
                        mutation = cov - ref

                        cov_data.append(cov)
                        cov_rows.append(i)
                        cov_cols.append(cell_idx)

                        mut_data.append(mutation / cov)
                        mut_rows.append(i)
                        mut_cols.append(cell_idx)
            except (TypeError, ValueError):
                continue

        if progress_interval and (i + 1) % progress_interval == 0:
            print(f"  Processed {i + 1:,} variants...")

    vcf.close()
    n_variants = len(row_ids)

    print(f"Parsed {n_variants:,} positions")
    print(f"Non-zero entries: {len(cov_data):,}")
    print(f"Sparsity: {100 * (1 - len(cov_data) / (n_variants * n_cells)):.2f}%")

    print("Building sparse matrices...")
    cov_sparse = sparse.csr_matrix(
        (cov_data, (cov_rows, cov_cols)),
        shape=(n_variants, n_cells), dtype=np.int32
    )
    mut_sparse = sparse.csr_matrix(
        (mut_data, (mut_rows, mut_cols)),
        shape=(n_variants, n_cells), dtype=np.float32
    )

    # Memory estimate
    cov_mem = (cov_sparse.data.nbytes + cov_sparse.indices.nbytes + cov_sparse.indptr.nbytes) / 1e9
    mut_mem = (mut_sparse.data.nbytes + mut_sparse.indices.nbytes + mut_sparse.indptr.nbytes) / 1e9
    dense_mem = n_variants * n_cells * 8 / 1e9  # What dense would use
    print(f"Memory usage: {cov_mem + mut_mem:.2f} GB (dense would be ~{dense_mem:.1f} GB)")

    print("Done!")
    print(f"Coverage matrix shape: {cov_sparse.shape}")
    print(f"Mutation rate matrix shape: {mut_sparse.shape}")

    return cov_sparse, mut_sparse, row_ids, samples


def read_cellsnp_base_vcf(vcf_file, include_oth=False, min_coverage=100):
    """
    Read cellSNP.base.vcf or cellSNP.base.vcf.gz file and generate a DataFrame 
    with gene, position, reference, coverage, and mutation rate information.
    
    Parameters:
    -----------
    vcf_file : str
        Path to cellSNP.base.vcf or cellSNP.base.vcf.gz file
        Automatically detects if file is gzipped or not
    include_oth : bool, default=False
        If True:  coverage = DP + OTH, mutrate = (AD + OTH) / (DP + OTH)
        If False: coverage = DP, mutrate = AD / DP
    min_coverage : int, default=100
        Minimum coverage threshold. Only positions with coverage >= min_coverage
        will be included in the final DataFrame.

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: ['gene', 'pos', 'ref', 'coverage', 'mutrate']
    """
    import pandas as pd
    import gzip
    import os

    print(f"Reading VCF file: {vcf_file}")
    print(f"Include OTH in calculations: {include_oth}")
    print(f"Minimum coverage threshold: {min_coverage}")

    # Check if file exists
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"File not found: {vcf_file}")
    
    # Determine if file is gzipped
    is_gzipped = vcf_file.endswith('.gz')
    
    # Lists to store data
    genes = []
    positions = []
    refs = []
    coverages = []
    mutrates = []
    
    # Counters
    n_records = 0
    n_filtered = 0
    
    # Open file appropriately based on compression
    if is_gzipped:
        print("Detected gzipped file")
        file_handle = gzip.open(vcf_file, 'rt')
    else:
        print("Detected uncompressed file")
        file_handle = open(vcf_file, 'r')
    
    try:
        for line in file_handle:
            # Skip header lines
            if line.startswith('#'):
                continue
            
            # Parse data line
            parts = line.strip().split('\t')
            
            if len(parts) < 8:
                continue
            
            gene = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            
            # Parse INFO field
            info_str = parts[7]
            info_dict = {}
            for item in info_str.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    try:
                        info_dict[key] = int(value)
                    except ValueError:
                        info_dict[key] = value
            
            # Extract values
            ad = info_dict.get('AD', 0)
            dp = info_dict.get('DP', 0)
            oth = info_dict.get('OTH', 0)
            
            # Calculate coverage based on include_oth
            if include_oth:
                coverage = dp + oth
            else:
                coverage = dp
            
            # Apply coverage threshold filter
            if coverage < min_coverage:
                n_filtered += 1
                continue
            
            # Calculate mutation rate
            if include_oth:
                mutrate = (ad + oth) / coverage if coverage > 0 else 0.0
            else:
                mutrate = ad / dp if dp > 0 else 0.0
            
            # Store values
            genes.append(gene)
            positions.append(pos)
            refs.append(ref)
            coverages.append(coverage)
            mutrates.append(mutrate)
            
            n_records += 1
    
    finally:
        # Always close the file handle
        file_handle.close()
    
    print(f"\nFiltering summary:")
    print(f"  Total records processed: {n_records + n_filtered}")
    print(f"  Records passing filter (coverage >= {min_coverage}): {n_records}")
    print(f"  Records filtered out: {n_filtered}")
    
    # Create DataFrame
    df = pd.DataFrame({
        'gene': genes,
        'pos': positions,
        'ref': refs,
        'coverage': coverages,
        'mutrate': mutrates
    })
    
    print("\nProcessing complete!")
    print(f"DataFrame shape: {df.shape}")
    
    if len(df) > 0:
        print(f"\nDataFrame summary:")
        print(f"  Genes: {df['gene'].nunique()} unique")
        print(f"  Positions: {len(df)}")
        print(f"  Coverage range: [{df['coverage'].min()}, {df['coverage'].max()}]")
        print(f"  Mean coverage: {df['coverage'].mean():.2f}")
        print(f"  Mutation rate range: [{df['mutrate'].min():.6f}, {df['mutrate'].max():.6f}]")
        print(f"  Mean mutation rate: {df['mutrate'].mean():.6f}")
    else:
        print("Warning: No records passed the coverage threshold!")
    
    return df


def read_window_vcf_to_sparse(vcf_file, progress_interval=100000):
    """
    Read window-based VCF file with per-cell DP:AD:MR format and generate sparse matrices.

    This function reads VCF files where each row represents a genomic window with
    per-cell coverage (DP), alternate counts (AD), and mutation rate (MR) information.

    Parameters:
    -----------
    vcf_file : str
        Path to VCF file (supports both .vcf and .vcf.gz)
        Format: gene, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT(DP:AD:MR), cell_data...
        Cell data format: "DP:AD:MR" (e.g., "103:0:0.0") or ".:.:." for missing
    progress_interval : int, default=100000
        Print progress every N variants. Set to 0 to disable progress output.

    Returns:
    --------
    tuple: (coverage_sparse, mutrate_sparse, position_info, col_names)
        - coverage_sparse: scipy.sparse.csr_matrix with DP (coverage) values
        - mutrate_sparse: scipy.sparse.csr_matrix with MR (mutation rate) values
        - position_info: pandas.DataFrame with columns ['gene', 'pos', 'ref', 'win_end']
        - col_names: list of cell/sample names
    """
    import gzip
    import pandas as pd
    import numpy as np
    from scipy import sparse
    import os

    print(f"Reading VCF file (sparse mode): {vcf_file}")

    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"File not found: {vcf_file}")

    is_gzipped = vcf_file.endswith('.gz')

    genes = []
    positions = []
    refs = []
    win_ends = []

    cov_data, cov_rows, cov_cols = [], [], []
    mut_data, mut_rows, mut_cols = [], [], []

    col_names = []
    n_cells = 0
    row_idx = 0

    if is_gzipped:
        file_handle = gzip.open(vcf_file, 'rt')
    else:
        file_handle = open(vcf_file, 'r')

    try:
        for line in file_handle:
            if line.startswith('##'):
                continue

            if line.startswith('#CHROM') or line.startswith('#'):
                parts = line.strip().split('\t')
                try:
                    format_idx = parts.index('FORMAT')
                    col_names = parts[format_idx + 1:]
                except ValueError:
                    col_names = parts[9:]
                n_cells = len(col_names)
                print(f"Found {n_cells} cells")
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            gene = parts[0]
            pos = int(parts[1])
            ref = parts[3]

            info_str = parts[7]
            win_end = None
            for item in info_str.split(';'):
                if item.startswith('WIN_END='):
                    try:
                        win_end = int(item.split('=')[1])
                    except ValueError:
                        pass
                    break

            genes.append(gene)
            positions.append(pos)
            refs.append(ref)
            win_ends.append(win_end)

            for cell_idx, cell_data in enumerate(parts[9:]):
                if cell_idx >= n_cells:
                    break

                if cell_data.startswith('.') or cell_data == '.:.:.':
                    continue

                cell_fields = cell_data.split(':')
                if len(cell_fields) < 3:
                    continue

                try:
                    dp_str, ad_str, mr_str = cell_fields[0], cell_fields[1], cell_fields[2]
                    if dp_str == '.' or mr_str == '.':
                        continue

                    dp_val = int(dp_str)
                    mr_val = float(mr_str)

                    if dp_val > 0:
                        cov_data.append(dp_val)
                        cov_rows.append(row_idx)
                        cov_cols.append(cell_idx)

                        mut_data.append(mr_val)
                        mut_rows.append(row_idx)
                        mut_cols.append(cell_idx)

                except (ValueError, IndexError):
                    continue

            row_idx += 1

            if progress_interval > 0 and row_idx % progress_interval == 0:
                print(f"  Processed {row_idx:,} variants...")

    finally:
        file_handle.close()

    n_variants = len(genes)
    print(f"Parsed {n_variants:,} positions")
    print(f"Non-zero entries: {len(cov_data):,}")

    if n_variants > 0 and n_cells > 0:
        print(f"Sparsity: {100 * (1 - len(cov_data) / (n_variants * n_cells)):.2f}%")

    cov_sparse = sparse.csr_matrix(
        (cov_data, (cov_rows, cov_cols)),
        shape=(n_variants, n_cells), dtype=np.int32
    )
    mut_sparse = sparse.csr_matrix(
        (mut_data, (mut_rows, mut_cols)),
        shape=(n_variants, n_cells), dtype=np.float32
    )

    position_info = pd.DataFrame({
        'gene': genes,
        'pos': positions,
        'ref': refs,
        'win_end': win_ends
    })
    position_info.index = [f"{g}.{p}" for g, p in zip(genes, positions)]

    print(f"Coverage matrix shape: {cov_sparse.shape}")
    print(f"Mutation rate matrix shape: {mut_sparse.shape}")

    return cov_sparse, mut_sparse, position_info, col_names

def read_paired_vcf_to_sparse(rRNA_win_path, steptwo_win_path, progress_interval=100000):
    """
    Read paired rRNA and steptwo VCF files and concatenate them.
    
    Parameters:
    -----------
    rRNA_win_path : str
        Path to rRNA VCF file
    steptwo_win_path : str
        Path to steptwo VCF file
    progress_interval : int, default=100000
        Print progress every N variants
        
    Returns:
    --------
    tuple: (coverage_sparse, mutrate_sparse, position_info, col_names)
        - coverage_sparse: scipy.sparse.csr_matrix with concatenated DP values
        - mutrate_sparse: scipy.sparse.csr_matrix with concatenated MR values
        - position_info: pandas.DataFrame with concatenated position info
        - col_names: list of cell/sample names (should be identical for both)
    """
    from scipy import sparse
    import pandas as pd
    
    print("=" * 80)
    print("Reading rRNA VCF file...")
    print("=" * 80)
    cov_rRNA, mut_rRNA, pos_rRNA, cells_rRNA = read_window_vcf_to_sparse(
        rRNA_win_path, progress_interval
    )
    
    print("\n" + "=" * 80)
    print("Reading steptwo VCF file...")
    print("=" * 80)
    cov_steptwo, mut_steptwo, pos_steptwo, cells_steptwo = read_window_vcf_to_sparse(
        steptwo_win_path, progress_interval
    )
    
    # Verify that cell names match
    if cells_rRNA != cells_steptwo:
        raise ValueError(
            f"Cell names don't match between rRNA and steptwo files!\n"
            f"rRNA cells: {len(cells_rRNA)}, steptwo cells: {len(cells_steptwo)}"
        )
    
    print("\n" + "=" * 80)
    print("Concatenating matrices...")
    print("=" * 80)
    
    # Concatenate sparse matrices vertically (row-wise)
    cov_combined = sparse.vstack([cov_rRNA, cov_steptwo], format='csr')
    mut_combined = sparse.vstack([mut_rRNA, mut_steptwo], format='csr')
    
    # Concatenate position info
    pos_combined = pd.concat([pos_rRNA, pos_steptwo], axis=0, ignore_index=False)
    
    print(f"Combined coverage matrix shape: {cov_combined.shape}")
    print(f"Combined mutation rate matrix shape: {mut_combined.shape}")
    print(f"Combined position info shape: {pos_combined.shape}")
    print(f"Total windows: {len(pos_combined):,}")
    print("=" * 80)
    
    return cov_combined, mut_combined, pos_combined, cells_rRNA


def read_multiple_paired_vcf_to_dict(sample_paths_dict, progress_interval=100000):
    """
    Read multiple samples with paired rRNA and steptwo VCF files.
    
    Parameters:
    -----------
    sample_paths_dict : dict
        Dictionary with sample names as keys and tuples of (rRNA_path, steptwo_path) as values
        Example: {
            'DM_PIP_DMSO_3in4': (rRNA_path, steptwo_path),
            'DM_PIP_NAIN3_3in4': (rRNA_path, steptwo_path)
        }
    progress_interval : int, default=100000
        Print progress every N variants
        
    Returns:
    --------
    dict: Dictionary with sample names as keys and tuples of 
          (coverage_sparse, mutrate_sparse, position_info, col_names) as values
    """
    results = {}
    
    for sample_name, (rRNA_path, steptwo_path) in sample_paths_dict.items():
        print("\n" + "=" * 80)
        print(f"PROCESSING SAMPLE: {sample_name}")
        print("=" * 80)
        
        try:
            cov, mut, pos, cells = read_paired_vcf_to_sparse(
                rRNA_path, steptwo_path, progress_interval
            )
            results[sample_name] = (cov, mut, pos, cells)
            
            print(f"\n✓ Successfully processed {sample_name}")
            print(f"  - Windows: {pos.shape[0]:,}")
            print(f"  - Cells: {len(cells):,}")
            print(f"  - Coverage non-zeros: {cov.nnz:,}")
            
        except Exception as e:
            print(f"\n✗ Error processing {sample_name}: {str(e)}")
            raise
    
    return results