import os
import re
import sys

def process_tex_file(filepath: str) -> str:
    """
    Recursively processes a LaTeX (.tex) file, replacing \input{} commands
    with the content of the specified input files.

    Args:
        filepath (str): The path to the LaTeX file to process.

    Returns:
        str: The fully processed content of the LaTeX file.
    """
    if not os.path.exists(filepath):
        print(f"Warning: File not found - {filepath}. Skipping input.", file=sys.stderr)
        return ""

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file {filepath}: {e}", file=os.stderr)
        return ""

    processed_content = []
    # Regex to find \input{filename} or \input{filename.tex}
    # It captures the filename inside the curly braces.
    input_pattern = re.compile(r'\\input\{(.*?)\}')

    for line in lines:
        match = input_pattern.search(line)
        if match:
            # Extract the filename from the match object
            input_filename = match.group(1)

            # Construct the full path for the input file
            # Assuming input files are relative to the current file's directory
            current_dir = os.path.dirname(filepath)
            
            # Add .tex extension if it's missing, as common in LaTeX
            if not input_filename.endswith('.tex'):
                input_filename += '.tex'
            
            full_input_path = os.path.join(current_dir, input_filename)

            print(f"Processing \\input: {full_input_path}")
            # Recursively process the input file and append its content
            processed_content.append(process_tex_file(full_input_path))
        else:
            # If no \input command, append the line as is
            processed_content.append(line)

    return "".join(processed_content)

if __name__ == "__main__":
    print("--- Starting LaTeX File Processing ---")
    output_content = process_tex_file("springer.tex")
    print("\n--- Processed Content ---")
    print(output_content)

    # Optionally, save the output to a new file
    output_filename = "springer_full.tex"
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(output_content)
    print(f"\nProcessed content saved to {output_filename}")
