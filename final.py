# Read the original LaTeX file
with open("springer.tex", "r", encoding="utf-8") as f:
    content = f.read()

# 1. Force insert the macro at the absolute top by replacing the first comment marker
macro_definition = "\\newcommand{\\mysecret}{2}\n"

# Replace the very first comment line with the macro + the comment line
content = content.replace("%Version 3.1 December 2024", f"{macro_definition}%Version 3.1 December 2024")

# 2. Anonymize the explicit names in the Funding/Declarations section
old_funding = (
    "\\item Funding Ondrej Šuch and Ali Haidar were  partially supported by grants "
    "VEGA 2/0172/22 Classification using ensembles of neural networks and "
    "VEGA 2/0056/25 Cycles and edge colorings of cubic graphs.\n"
    "Peter Novotný did not receive support from any organization for the submitted work."
)
new_funding = "\\item Funding: [Anonymized for double-blind review]."

content = content.replace(old_funding, new_funding)

# Write the fully fixed version directly to sp2.tex
with open("sp2.tex", "w", encoding="utf-8") as f:
    f.write(content)

print("Success! File saved as 'sp2.tex'.")
