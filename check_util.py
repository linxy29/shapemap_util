def read_dot_bracket(file_path, gene_name):
    """
    Read a file containing dot-bracket notation and extract the structure for a specific gene.
    
    Parameters:
    -----------
    file_path : str
        Path to the file containing dot-bracket notation (e.g., 'references/Dot bracket for ribo for mouse.txt')
    gene_name : str
        Name of the gene to extract (e.g., '18S', '12S', '16S')
    
    Returns:
    --------
    str
        Dot-bracket notation string for the specified gene, or None if gene not found
    
    Example:
    --------
    >>> structure = read_dot_bracket('references/Dot bracket for ribo for mouse.txt', '18S')
    >>> print(structure[:50])  # Print first 50 characters
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Look for the gene name (it will be on a line starting with '>')
        target_header = f'>{gene_name}'
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Check if this line contains the gene name
            if line.startswith(target_header):
                # The dot-bracket notation should be on the next line(s)
                dot_bracket = ''
                
                # Read subsequent lines until we hit another header or end of file
                for j in range(i + 1, len(lines)):
                    next_line = lines[j].strip()
                    
                    # Stop if we encounter another gene header
                    if next_line.startswith('>'):
                        break
                    
                    # Accumulate the dot-bracket notation
                    dot_bracket += next_line
                
                return dot_bracket if dot_bracket else None
        
        # Gene not found
        print(f"Gene '{gene_name}' not found in file.")
        return None
    
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None