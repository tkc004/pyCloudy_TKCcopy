from pathlib import Path


def test_docs_pages_and_example_exports_exist():
    root = Path(__file__).resolve().parents[1]

    required_docs = [
        root / "docs" / "README.md",
        root / "docs" / "API.md",
        root / "docs" / "CloudyInput.md",
        root / "docs" / "index.html",
        root / "python_examples" / "README.md",
    ]

    for path in required_docs:
        assert path.exists(), f"missing documentation file: {path}"


def test_notebooks_have_python_exports():
    root = Path(__file__).resolve().parents[1]
    notebook_dir = root / "pyCloudy" / "docs"
    export_dir = root / "python_examples"

    notebook_stems = sorted(p.stem for p in notebook_dir.glob("*.ipynb"))

    assert notebook_stems, "expected at least one notebook in pyCloudy/docs"
    for stem in notebook_stems:
        export = export_dir / stem / f"{stem}.py"
        assert export.exists(), (
            f"python_examples should contain {stem}/{stem}.py for notebook {stem}.ipynb"
        )


def test_quickstart_script_exists():
    root = Path(__file__).resolve().parents[1]
    quickstart = root / "python_examples" / "quickstart" / "quickstart.py"

    assert quickstart.exists(), "missing quickstart example script"
