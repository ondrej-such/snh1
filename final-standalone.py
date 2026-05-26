import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "springer.tex"
OUTPUT = ROOT / "sp-standalone.tex"

INPUT_PATTERN = re.compile(r"\\input\{([^{}]+\.tex)\}")
IFTHEN_PATTERN = re.compile(r"\\ifthenelse\{\\not\\equal\{\\mysecret\}\{(\d)\}\}")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def is_commented_command(text: str, pos: int) -> bool:
    line_start = text.rfind("\n", 0, pos) + 1
    prefix = text[line_start:pos]
    return prefix.lstrip().startswith("%")


def find_matching_brace(text: str, open_index: int) -> int:
    if open_index >= len(text) or text[open_index] != "{":
        raise ValueError(f"Expected '{{' at position {open_index}")

    depth = 0
    i = open_index
    while i < len(text):
        ch = text[i]

        if ch == "%":
            newline = text.find("\n", i)
            if newline == -1:
                return len(text) - 1
            i = newline + 1
            continue

        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i

        i += 1

    raise ValueError("Unmatched brace in LaTeX content")


def expand_inputs(text: str, base_dir: Path) -> str:
    result = []
    i = 0

    while i < len(text):
        match = INPUT_PATTERN.search(text, i)
        if not match:
            result.append(text[i:])
            break

        start, end = match.span()
        result.append(text[i:start])

        if is_commented_command(text, start):
            result.append(text[start:end])
        else:
            rel_path = match.group(1)
            file_path = (base_dir / rel_path).resolve()
            result.append(read_text(file_path))

        i = end

    expanded = "".join(result)
    if expanded == text:
        return expanded
    return expand_inputs(expanded, base_dir)


def process_ifthenelse(text: str) -> str:
    result = []
    i = 0

    while i < len(text):
        match = IFTHEN_PATTERN.search(text, i)
        if not match:
            result.append(text[i:])
            break

        start = match.start()
        end = match.end()

        result.append(text[i:start])

        if is_commented_command(text, start):
            result.append(text[start:end])
            i = end
            continue

        secret_value = match.group(1)

        if end >= len(text) or text[end] != "{":
            raise ValueError(f"Missing true branch near position {start}")

        true_open = end
        true_close = find_matching_brace(text, true_open)
        true_content = text[true_open + 1:true_close]

        cursor = true_close + 1

        if cursor < len(text) and text[cursor] == "{":
            false_close = find_matching_brace(text, cursor)
            cursor = false_close + 1

        if secret_value != "2":
            result.append(true_content)

        i = cursor

    return "".join(result)


def main() -> None:
    content = read_text(SOURCE)
    content = expand_inputs(content, ROOT)
    content = process_ifthenelse(content)
    write_text(OUTPUT, content)
    print(f"Success! File saved as '{OUTPUT.name}'.")


if __name__ == "__main__":
    main()