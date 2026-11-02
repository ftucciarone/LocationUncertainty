#!/usr/bin/env python3
"""
Write docs/py-kits.html, a reference with one page per kit, module and function, laid out like
the MathWorks function pages: syntax, description and formulas, examples, input and output
arguments, NEMO implementation, see also. Everything comes from the docstrings, the examples
from tests/test_all_routines.ipynb, and the NEMO lines from the source trees in paths.toml
(nemo_src, newest first), linked to forge.nemo-ocean.eu at the release that has them.
The page holds no local path, so it can be published as it is.

    python docs/make_html.py                # from py-kits; or name kits, e.g. orcakit datakit
"""
import ast
import html
import json
import os
import re
import subprocess
import sys
import textwrap

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from orcakit import paths  # noqa: E402

FORGE = "https://forge.nemo-ocean.eu/nemo/nemo/-/blob"
KATEX = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9"
KITS = ("datakit", "orcakit", "noisekit", "modredkit", "diagkit", "gyrekit")
SECTIONS = ("Parameters:", "Returns:", "Formula:", "NEMO:", "See also:")
STYLE = """
:root { --bg: #ffffff; --fg: #1d1d1f; --muted: #66666c; --line: #dcdce0; --code: #f3f4f6; --accent: #1f5fa8; --warn: #a3470d; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { --bg: #17181a; --fg: #e8e8ea; --muted: #a0a0a8; --line: #33353a; --code: #212327; --accent: #7fb0ea; --warn: #e7a06a; }
}
body { margin: 0; background: var(--bg); color: var(--fg); font: 15px/1.55 system-ui, sans-serif; display: flex; }
nav { width: 230px; flex: none; border-right: 1px solid var(--line); padding: 16px; position: sticky; top: 0; height: 100vh; box-sizing: border-box; overflow-y: auto; font-size: 14px; }
nav a { display: block; color: var(--fg); text-decoration: none; padding: 2px 0 2px 12px; }
nav summary { cursor: pointer; padding: 3px 0; font-weight: 600; }
nav a:hover, a { color: var(--accent); }
main { flex: 1; min-width: 0; max-width: 960px; padding: 16px 32px 64px; }
section { display: none; }
section:target { display: block; }
main:not(:has(section:target)) #home { display: block; }
.crumb { color: var(--muted); font-size: 13px; }
h1 { font-size: 26px; margin: 4px 0; font-family: ui-monospace, monospace; }
h2 { font-size: 18px; margin-top: 32px; border-bottom: 1px solid var(--line); padding-bottom: 4px; }
.summary { font-size: 16px; margin-top: 0; }
.toc { font-size: 13px; color: var(--muted); }
.toc a { margin-right: 12px; }
pre, code { font-family: ui-monospace, monospace; font-size: 13px; }
pre { background: var(--code); padding: 10px 12px; border-radius: 6px; overflow-x: auto; }
.syntax code { display: block; background: var(--code); padding: 6px 10px; border-radius: 6px; margin: 6px 0; }
.math { overflow-x: auto; }
dl.args dt { font-family: ui-monospace, monospace; font-weight: 600; margin-top: 14px; }
dl.args dt span { font-family: system-ui, sans-serif; font-weight: 400; color: var(--muted); }
dl.args dd { margin: 4px 0 0 18px; }
.nemo a { font-family: ui-monospace, monospace; font-size: 13px; }
.missing { color: var(--warn); }
table { border-collapse: collapse; width: 100%; }
td { border-top: 1px solid var(--line); padding: 6px 8px; vertical-align: top; }
td:first-child { font-family: ui-monospace, monospace; white-space: nowrap; }
@media (max-width: 760px) { body { display: block; } nav { width: auto; height: auto; position: static; border-right: 0; border-bottom: 1px solid var(--line); } main { padding: 16px; } }
@media print { nav { display: none; } section { display: block !important; break-before: page; } }
"""


def parse_doc(doc, name):
    r"""
    parse_doc: split a docstring in the py-kits style into its parts.

        Parameters:
            doc  (str): Docstring, dedented
            name (str): Function name, stripped from the summary

        Returns:
            summary   (str): First line, without 'name: '
            prose    (list): Paragraphs of the description, as lists of lines
            sections (dict): {'Parameters:': lines, ...}, lines dedented by 4

        Formula:
            \text{doc} = \text{name: summary} \oplus \text{prose} \oplus \{\text{Parameters, Returns, Formula, NEMO, See also}\}

        See also:
            entries
    """
    lines = (doc or "").splitlines()
    first = 1
    while first < len(lines) and lines[first].strip() and not lines[first].startswith("    "):
        first += 1
    head = " ".join(l.strip() for l in lines[:first])
    summary = head[len(name) + 2:] if head.startswith(name + ":") else head
    prose, sections, current = [[]], {}, None
    for line in lines[first:]:
        if line.strip() in SECTIONS and line.startswith("    ") and not line.startswith("        "):
            current = line.strip()
            sections[current] = []
        elif current:
            sections[current].append(line[4:])
        elif line.strip():
            prose[-1].append(line)
        elif prose[-1]:
            prose.append([])
    return summary, [p for p in prose if p], sections


def entries(lines):
    r"""
    entries: the 'name (type): description' entries of a Parameters or Returns section.

        Parameters:
            lines (list): Section lines, entries indented by 4, continuations further

        Returns:
            out (list): (name, type, description) per entry; type is '' when the entry has none

        Formula:
            \text{line} \sim \texttt{name (type): description},\ \text{continuations appended}

        See also:
            parse_doc
    """
    out = []
    for line in lines:
        m = re.match(r"^    (\S.*?)\s*\(([^()]*)\):\s*(.*)$", line)
        if m:
            out.append([m.group(1), m.group(2), m.group(3)])
        elif line.strip() and out:
            out[-1][2] += " " + line.strip()
        elif line.strip():
            out.append(["", "", line.strip()])
    return out


def defaults(node):
    r"""
    defaults: the default value of each argument that has one, as written in the source.

        Parameters:
            node (FunctionDef): Function

        Returns:
            out (dict): {argument: source text of its default}

        Formula:
            \{a : \text{source text of its default}\}

        See also:
            function_page
    """
    a = node.args
    pos = a.args[len(a.args) - len(a.defaults):]
    out = {arg.arg: ast.unparse(d) for arg, d in zip(pos, a.defaults)}
    out.update({arg.arg: ast.unparse(d) for arg, d in zip(a.kwonlyargs, a.kw_defaults) if d is not None})
    return out


def examples(notebook):
    r"""
    examples: the test code and results of the test notebook, by routine.

    Each test block is matched with the header it prints, which names every routine it tests, even names
    built in a loop.

        Parameters:
            notebook (str): Path of the executed notebook

        Returns:
            out (dict): {'kit.module.function': (code, result lines)}

        Formula:
            \text{block}_n \leftrightarrow \text{header}_n \text{ in the cell's output},\qquad \text{routine} \mapsto (\text{code}, \text{results})

        See also:
            function_page
    """
    out = {}
    with open(notebook) as f:
        cells = json.load(f)["cells"]
    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"]).splitlines()
        text = "".join(o.get("text", "") if isinstance(o.get("text", ""), str) else "".join(o["text"])
                       for o in cell.get("outputs", []))
        groups = []
        for block in re.split(r"^== ", text, flags=re.M)[1:]:
            head, _, body = block.partition("\n")
            lines = [l.strip() for l in body.splitlines() if l.strip().startswith(("PASS", "FAIL", "ERROR"))]
            groups.append(([n.split(" (")[0] for n in head.split(", ")], lines))
        blocks = []
        for i, line in enumerate(src):
            m = re.match(r"^(\s*)with test\((.*)\):", line)
            if not m:
                continue
            body = [line]
            for nxt in src[i + 1:]:
                if nxt.strip() and len(nxt) - len(nxt.lstrip()) <= len(m.group(1)):
                    break
                body.append(nxt)
            blocks.append("\n".join(body).rstrip())
        # The printed headers name every routine a block tests, even names built in a loop
        for code, (names, lines) in zip(blocks, groups):
            for name in names:
                out[name] = (code, lines)
    return out


class Usage:
    """Where each routine of a kit is called or passed in the library: the kits and preprocessing/."""

    def __init__(self, kits, roots):
        r"""
        __init__: index every use of the kits' routines in the library.

        A use is a call or a reference, such as exchange=orca.halo, through any import alias, or by bare name
        inside the routine's own module. Tests are left out.

            Parameters:
                kits  (list): Kits whose routines are looked for, e.g. ['orcakit', 'datakit']
                roots (list): (label, folder) pairs of Python code to search

            Returns:
                None

            Formula:
                \text{uses}[\text{routine}][(\text{code base}, \text{file}, \text{caller})] = (\text{caller's summary}, \text{snippets})

            See also:
                Usage.index
        """
        self.kits, self.uses = tuple(k + "." for k in kits), {}
        for label, root in roots:
            for folder, _, names in os.walk(root):
                if "tests" in folder.split(os.sep) or "__pycache__" in folder:
                    continue
                for f in sorted(names):
                    if f.endswith(".py"):
                        self.index(label, root, os.path.join(folder, f))

    def index(self, label, root, path):
        r"""
        index: record the uses in one file.

            Parameters:
                label (str): Name shown for the code base, e.g. 'py-kits'
                root  (str): Its folder
                path  (str): Python file

            Returns:
                None

            Formula:
                \text{for each name or attribute resolving to a kit routine, inside a function: record the statement}

            See also:
                Usage.__init__
        """
        with open(path) as f:
            src = f.read()
        tree = ast.parse(src)
        rel = os.path.relpath(path, root)
        here = rel[:-3].replace(os.sep, ".")
        alias = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                module = node.module if node.level == 0 else ".".join(here.split(".")[:-node.level] + [node.module])
                for a in node.names:
                    alias[a.asname or a.name] = f"{module}.{a.name}"
            elif isinstance(node, ast.Import):
                for a in node.names:
                    alias[a.asname or a.name.split(".")[0]] = a.name if a.asname else a.name.split(".")[0]

        def target(node):
            r"""
            target: the full name a name or attribute refers to, through the file's imports.

                Parameters:
                    node (AST): Name or Attribute node

                Returns:
                    name (str): Dotted full name, e.g. 'orcakit.orca.halo', or None

                Formula:
                    \text{alias}[\text{base}] \oplus \text{attributes},\ \text{or the module's own name for a bare name}

                See also:
                    Usage.index
            """
            if isinstance(node, ast.Attribute):
                base, parts = node.value, [node.attr]
                while isinstance(base, ast.Attribute):
                    parts.insert(0, base.attr)
                    base = base.value
                if isinstance(base, ast.Name) and base.id in alias:
                    return ".".join([alias[base.id]] + parts)
            elif isinstance(node, ast.Name):
                if node.id in alias:
                    return alias[node.id]
                if here.startswith(self.kits):
                    return f"{here}.{node.id}"
            return None

        def visit(node, stack, stmt):
            r"""
            visit: walk the tree, tracking the enclosing functions and the current statement.

                Parameters:
                    node  (AST): Node whose children are visited
                    stack (list): Enclosing function and class nodes, outermost first
                    stmt  (AST): Statement the node belongs to

                Returns:
                    None; uses are recorded in self.uses

                Formula:
                    \text{depth-first over the syntax tree}

                See also:
                    Usage.index
            """
            for c in ast.iter_child_nodes(node):
                s = c if isinstance(c, ast.stmt) else stmt
                if isinstance(c, (ast.FunctionDef, ast.ClassDef)):
                    visit(c, stack + [c], s)
                    continue
                if isinstance(c, (ast.Attribute, ast.Name)) and not isinstance(getattr(c, "ctx", None), ast.Store):
                    t = target(c)
                    funcs = [n for n in stack if isinstance(n, ast.FunctionDef)]
                    if t and t.startswith(self.kits) and funcs and t != f"{here}.{'.'.join(n.name for n in stack)}":
                        caller = next((n for n in reversed(funcs) if ast.get_docstring(n)), funcs[-1])
                        qual = ".".join(n.name for n in stack)
                        lines = (ast.get_source_segment(src, s) or "").splitlines()
                        snippet = "\n".join(lines[:4]) + ("\n..." if len(lines) > 4 else "")
                        entry = self.uses.setdefault(t, {}).setdefault((label, rel, qual), [parse_doc(ast.get_docstring(caller), caller.name)[0], []])
                        if snippet not in entry[1]:
                            entry[1].append(snippet)
                visit(c, stack, s)
        visit(tree, [], None)


class Nemo:
    """Find quoted lines in NEMO source trees, newest release first."""

    def __init__(self, roots):
        r"""
        __init__: the NEMO source trees and their release tags.

        Trees that are not git clones at a release tag are skipped with a warning, since their lines could not
        be linked to the forge.

            Parameters:
                roots (list): Folders of NEMO git clones; their order does not matter

            Returns:
                None

            Formula:
                \text{trees sorted by release, newest first}

            See also:
                Nemo.find
        """
        self.trees = []
        for root in roots:
            try:
                tag = subprocess.run(["git", "-C", root, "describe", "--tags", "--exact-match"],
                                     capture_output=True, text=True, check=True).stdout.strip()
            except (OSError, subprocess.CalledProcessError):
                print(f"warning: {root} is not a NEMO clone at a release tag; skipped", file=sys.stderr)
                continue
            self.trees.append((tuple(int(x) for x in re.findall(r"\d+", tag)), tag, root))
        self.trees.sort(reverse=True)
        self.missing = []

    def find(self, path, quote):
        r"""
        find: the newest release whose file holds a quoted line.

        Citing the newest release that still has the line keeps links current, and still finds routines that
        were removed in a later release.

            Parameters:
                path  (str): File, relative to the NEMO root
                quote (str): Text of the line, compared with spaces collapsed

            Returns:
                ref (tuple): (tag, line number, url), or None if no tree has it

            Formula:
                \min_{\text{release, newest first}} \{\, n : \text{quote} \subseteq \text{line}_n \,\}

            See also:
                function_page
        """
        want = " ".join(quote.split())
        for _, tag, root in self.trees:
            try:
                with open(os.path.join(root, path), errors="replace") as f:
                    for n, line in enumerate(f, 1):
                        if want in " ".join(line.split()):
                            return tag, n, f"{FORGE}/{tag}/{path}#L{n}"
            except OSError:
                continue
        self.missing.append(f"{path}: {quote}")
        return None


def capital(text):
    r"""
    capital: text with its first letter in upper case, the rest unchanged.

        Parameters:
            text (str): Text

        Returns:
            out (str): Text, first letter capitalised

        Formula:
            t_0 \to \mathrm{upper}(t_0)

        See also:
            function_page
    """
    return text[:1].upper() + text[1:]


def clean_doc(doc):
    r"""
    clean_doc: a module or package docstring without its title lines.

    Old docstrings start with a title such as 'Module orca', 'mpi_orca.py' and a line of '=', a path,
    or a creation date and author; the page already shows the name, so those are dropped.

        Parameters:
            doc (str): Docstring, or None

        Returns:
            text (str): Docstring from its first real paragraph

        Formula:
            \text{drop leading lines matching a title, an underline, a path, a date or an author}

        See also:
            render_text
    """
    lines = (doc or "").splitlines()
    title = re.compile(r"^(Module \w+|[\w/.-]+\.py|=+|-+|Created on .*|@author.*|)$")
    while lines and title.match(lines[0].strip()):
        lines.pop(0)
    return "\n".join(lines)


def render_text(text, base=None, known=()):
    r"""
    render_text: HTML of a free-form docstring: paragraphs, indented blocks, and a Modules list.

    A 'Modules:' section, with one 'name  description' entry per line, becomes a table whose names link
    to the module pages.

        Parameters:
            text   (str): Docstring, cleaned
            base   (str): Page id the module names are relative to, e.g. 'orcakit'
            known  (set): Page ids that exist, for the links

        Returns:
            html (str): HTML

        Formula:
            \text{paragraph} \to \langle p\rangle,\quad \text{indented block} \to \langle pre\rangle,\quad \text{Modules} \to \langle table\rangle

        See also:
            clean_doc
    """
    out, para, block, rows, in_modules = [], [], [], [], False

    def flush():
        if para:
            out.append(f"<p>{esc(' '.join(l.strip() for l in para))}</p>")
            para.clear()
        if block:
            out.append(f"<pre>{esc(textwrap.dedent(chr(10).join(block)).strip(chr(10)))}</pre>")
            block.clear()

    for line in text.splitlines() + [""]:
        if line.strip() == "Modules:":
            flush()
            in_modules = True
            continue
        if in_modules:
            m = re.match(r"^\s{2,}(\S+)\s{2,}(.*)$", line)
            if m:
                rows.append([m.group(1), m.group(2)])
                continue
            if line.startswith("      ") and line.strip() and rows:
                rows[-1][1] += " " + line.strip()
                continue
            if line.strip():
                in_modules = False
            else:
                continue
        if not line.strip():
            if block:
                block.append("")
            else:
                flush()
        elif line.startswith("    "):
            if para:
                flush()
            block.append(line)
        else:
            if block:
                flush()
            para.append(line)
    flush()
    if rows:
        cells = []
        for name, desc in rows:
            target = f"{base}.{name}" if base else name
            label = f'<a href="#{target}">{esc(name)}</a>' if target in known else esc(name)
            cells.append(f"<tr><td>{label}</td><td>{esc(desc)}</td></tr>")
        out.append("<h2>Modules</h2><table>" + "".join(cells) + "</table>")
    return "\n".join(out)


def first_sentence(doc):
    r"""
    first_sentence: the opening line of a cleaned docstring, for tables.

        Parameters:
            doc (str): Docstring, or None

        Returns:
            line (str): Its first paragraph, joined on one line

        Formula:
            \text{first paragraph of clean\_doc(doc)}

        See also:
            clean_doc
    """
    text = clean_doc(doc).strip().split("\n\n")[0]
    return " ".join(l.strip() for l in text.splitlines())


def esc(text):
    r"""
    esc: text made safe to place in HTML.

        Parameters:
            text (str): Text

        Returns:
            safe (str): Text with &, < and > escaped

        Formula:
            \& \to \&amp;,\quad < \ \to \&lt;,\quad > \ \to \&gt;

        See also:
            function_page
    """
    return html.escape(text, quote=False)


def link(name, kit, module, known):
    r"""
    link: an HTML link for a name in a 'See also' list.

    The name is looked up as written, then within the page's kit and module, then as the end of any
    page id of the page's kit, then of any kit.

        Parameters:
            name    (str): 'function', 'module.function' or 'kit.module.function'
            kit     (str): Kit of the page the link is on
            module  (str): Module of that page, within its kit
            known   (set): Page ids 'kit.module.function' that exist

        Returns:
            a (str): Link to the page, or the name in code style if no page matches

        Formula:
            \text{first existing of } n,\ k.n,\ k.m.n,\ \text{the unique id of } k \text{ ending in } .n,\ \text{the unique id ending in } .n

        See also:
            function_page
    """
    target = next((t for t in (name, f"{kit}.{name}", f"{kit}.{module}.{name}") if t in known), None)
    if target is None:
        head, _, rest = name.partition(".")
        scopes = [(f"{head}.", rest)] if head in KITS and rest else []
        for scope, tail in scopes + [(f"{kit}.", name), ("", name)]:
            ends = [k for k in known if k.startswith(scope) and k.endswith("." + tail)]
            if len(ends) == 1:
                target = ends[0]
                break
    return f'<a href="#{target}"><code>{esc(name)}</code></a>' if target else f"<code>{esc(name)}</code>"


def function_page(kit, module, qual, node, nemo, tests, known, usage):
    r"""
    function_page: the HTML section of one function.

    Summary, description and formula, syntax, arguments, example, NEMO lines, uses in the library and
    related routines, in that order.

        Parameters:
            kit      (str): Kit name
            module   (str): Module name
            qual     (str): Qualified name within the module, e.g. 'Decomp.slices'
            node (FunctionDef): Function
            nemo    (Nemo): NEMO source finder
            tests   (dict): From examples
            usage  (Usage): Where the kit's routines are used
            known    (set): Page ids 'kit.module.function' that exist, for links

        Returns:
            text (str): HTML

        Formula:
            \text{section}(\text{id} = \text{kit.module.qual})

        See also:
            parse_doc, examples, Usage, Nemo
    """
    name = node.name
    summary, prose, sec = parse_doc(ast.get_docstring(node), name)
    args = [a.arg for a in node.args.args + node.args.kwonlyargs if a.arg not in ("self", "cls")]
    if node.args.vararg:
        args.append("*" + node.args.vararg.arg)
    dflt = defaults(node)
    outs = [e[0] for e in entries(sec.get("Returns:", [])) if e[0] and e[0] != "None"]
    lhs = (", ".join(outs) + " = ") if outs else ""
    call = qual if "." not in qual or qual.split(".")[0][:1].isupper() else name
    required = [a for a in args if a not in dflt]
    forms = [f"{lhs}{call}({', '.join(required)})"]
    if len(required) < len(args):
        forms.append(f"{lhs}{call}({', '.join(required + [f'{a}={dflt[a]}' for a in args if a in dflt])})")
    full = pid = f"{kit}.{module}.{qual}"
    parts = [f'<section id="{pid}">',
             f'<div class="crumb"><a href="#home">py-kits</a> &rsaquo; <a href="#{kit}">{kit}</a> &rsaquo; <a href="#{kit}.{module}">{module}</a> &rsaquo; {esc(qual)}</div>',
             f"<h1>{esc(qual)}</h1>",
             f'<p class="summary">{esc(capital(summary))}</p>',
             "<h2>Description</h2>"]
    parts += [f"<p>{esc(' '.join(l.strip() for l in p))}</p>" if not p[0].startswith("    ")
              else f"<pre>{esc(chr(10).join(p))}</pre>" for p in prose]
    formula = [l.strip() for l in sec.get("Formula:", []) if l.strip()]
    parts += [f'<div class="math">$${esc(f)}$$</div>' for f in formula]
    parts += ["<h2>Syntax</h2>", '<div class="syntax">' + "".join(f"<code>{esc(f)}</code>" for f in forms) + "</div>"]
    for title, key in (("Input Arguments", "Parameters:"), ("Output Arguments", "Returns:")):
        ents = [e for e in entries(sec.get(key, [])) if e[0]]
        parts.append(f"<h2>{title}</h2>")
        if not ents:
            parts.append("<p>None.</p>")
            continue
        parts.append('<dl class="args">')
        for arg, typ, desc in ents:
            extra = f" &mdash; default <code>{esc(dflt[arg])}</code>" if arg in dflt else ""
            parts.append(f"<dt>{esc(arg)} <span>&mdash; {esc(typ)}{extra}</span></dt><dd>{esc(desc)}</dd>")
        parts.append("</dl>")
    parts.append("<h2>Example</h2>")
    if full in tests:
        code, results = tests[full]
        calls = [l.strip() for l in code.splitlines()[1:] if re.search(rf"\b{re.escape(name)}\(", l)]
        if calls:
            parts.append("<pre>" + esc("\n".join(calls[:12])) + "</pre>")
        parts.append(f"<details{'' if calls else ' open'}><summary>Whole test</summary><pre>{esc(code)}</pre></details>")
        if results:
            parts.append("<p>Result in <code>tests/test_all_routines.ipynb</code>:</p><pre>" + esc("\n".join(results)) + "</pre>")
    else:
        parts.append('<p class="toc">No test of its own; see the routines that call it.</p>')
    refs = [l.strip() for l in sec.get("NEMO:", []) if l.strip()]
    if refs:
        parts.append('<h2>NEMO implementation</h2><div class="nemo">')
        for ref in refs:
            path, _, quote = ref.partition(": ")
            found = nemo.find(path, quote)
            if found:
                tag, n, url = found
                parts.append(f'<p><a href="{url}">{esc(path)}, line {n}</a> (NEMO {tag})</p><pre>{esc(quote)}</pre>')
            else:
                parts.append(f'<p class="missing">{esc(path)}: line not found in the NEMO releases searched</p><pre>{esc(quote)}</pre>')
        parts.append("</div>")
    users = usage.uses.get(full, {})
    if users:
        parts.append("<h2>Used in the library by</h2>")
        for (label, rel, caller), (why, snippets) in sorted(users.items()):
            what = capital(why) if why else "No description."
            parts.append(f"<p><code>{esc(caller)}</code> in <code>{esc(label)}/{esc(rel)}</code>: {esc(what)}</p>")
            parts.append("<pre>" + esc("\n\n".join(textwrap.dedent(x) for x in snippets[:3])) + "</pre>")
    see = [s.strip() for l in sec.get("See also:", []) for s in l.split(",") if s.strip()]
    if see:
        parts.append("<h2>See also</h2><p>" + " | ".join(link(s, kit, module, known) for s in see) + "</p>")
    return "\n".join(parts + ["</section>"])


def functions(tree):
    r"""
    functions: every function of a module, nested ones included, with their qualified names.

        Parameters:
            tree (Module): Parsed module

        Returns:
            out (list): (qualified name, FunctionDef) pairs, in file order

        Formula:
            \{(\text{qual}, \text{node})\} \text{ in file order}

        See also:
            main
    """
    out = []

    def visit(node, prefix):
        r"""
        visit: walk the tree, collecting functions with their qualified names.

            Parameters:
                node   (AST): Node whose children are visited
                prefix (str): Qualified name of the enclosing classes and functions, with a trailing dot

            Returns:
                None; functions are appended to out

            Formula:
                \text{depth-first over the syntax tree}

            See also:
                functions
        """
        for c in ast.iter_child_nodes(node):
            if isinstance(c, (ast.FunctionDef, ast.ClassDef)):
                if isinstance(c, ast.FunctionDef):
                    out.append((prefix + c.name, c))
                visit(c, prefix + c.name + ".")
            else:
                visit(c, prefix)
    visit(tree, "")
    return out


def main():
    """
    main: write docs/py-kits.html, the reference of every kit; run from py-kits.

        Parameters:
            None; kits can be listed on the command line, all by default

        Returns:
            None; exits with status 1 if a NEMO quote was not found

        Formula:
            \text{one page per kit, module and function}

        See also:
            function_page, Usage, Nemo, examples
    """
    kits = sys.argv[1:] or list(KITS)
    nemo = Nemo([r for r in paths.get("nemo_src").split(":") if os.path.isdir(r)])
    tests = examples(os.path.join("tests", "test_all_routines.ipynb"))
    usage = Usage(kits, [("py-kits", "."),
                         ("preprocessing", os.path.join("..", "preprocessing"))])
    mods, packages = {}, {}
    for kit in kits:
        for folder, _, names in sorted(os.walk(kit)):
            if "__pycache__" in folder:
                continue
            for f in sorted(names):
                if not f.endswith(".py"):
                    continue
                with open(os.path.join(folder, f)) as src:
                    tree = ast.parse(src.read())
                if f == "__init__.py":
                    sub = os.path.relpath(folder, kit).replace(os.sep, ".")
                    packages[kit if sub == "." else f"{kit}.{sub}"] = ast.get_docstring(tree)
                else:
                    mods[(kit, os.path.relpath(os.path.join(folder, f[:-3]), kit).replace(os.sep, "."))] = tree
    known = {f"{k}.{m}.{q}" for (k, m), t in mods.items() for q, _ in functions(t)}
    known |= {f"{k}.{m}" for k, m in mods} | set(packages)
    nav, pages, count = [], [], 0
    for kit in kits:
        kmods = [(m, t) for (k, m), t in mods.items() if k == kit]
        nav.append(f'<details><summary><a href="#{kit}" style="display:inline;padding:0">{kit}</a></summary>'
                   + "".join(f'<a href="#{kit}.{m}">{m}</a>' for m, _ in kmods) + "</details>")
        for pkg in sorted(p for p in packages if p == kit or p.startswith(kit + ".")):
            crumb = '<a href="#home">py-kits</a>' + ("" if pkg == kit else f' &rsaquo; <a href="#{kit}">{kit}</a>')
            pages.append(f'<section id="{pkg}"><div class="crumb">{crumb}</div><h1>{pkg}</h1>'
                         f'<p class="summary">{esc(first_sentence(packages[pkg]))}</p>'
                         + render_text(clean_doc(packages[pkg]).split("\n\n", 1)[-1] if "\n\n" in clean_doc(packages[pkg]) else "", pkg, known)
                         + "</section>")
        for m, t in kmods:
            funcs = functions(t)
            rows = "".join(f'<tr><td><a href="#{kit}.{m}.{q}">{esc(q)}</a></td><td>{esc(capital(parse_doc(ast.get_docstring(n), n.name)[0]))}</td></tr>'
                           for q, n in funcs)
            doc = clean_doc(ast.get_docstring(t))
            rest = doc.split("\n\n", 1)[1] if "\n\n" in doc else ""
            parent = f"{kit}.{m.rsplit('.', 1)[0]}" if "." in m else kit
            pages.append(f'<section id="{kit}.{m}"><div class="crumb"><a href="#home">py-kits</a> &rsaquo; <a href="#{kit}">{kit}</a>'
                         + (f' &rsaquo; <a href="#{parent}">{parent.split(".", 1)[1]}</a>' if parent != kit else "")
                         + f'</div><h1>{m}</h1><p class="summary">{esc(first_sentence(doc))}</p>{render_text(rest)}'
                         f'<h2>Functions</h2><table>{rows}</table></section>')
            pages += [function_page(kit, m, q, n, nemo, tests, known, usage) for q, n in funcs]
            count += len(funcs)
    releases = ", ".join(t for _, t, _ in nemo.trees)
    with open(os.path.join("tests", "test_all_routines.ipynb")) as f:
        totals = re.findall(r"(\d+) test groups, (\d+) checks", f.read())
    untested = sorted(k for k in known if k.count(".") >= 2 and k not in tests and k not in packages
                      and not any(k == f"{kk}.{mm}" for kk, mm in mods))
    home = ('<section id="home"><h1>py-kits</h1><p class="summary">Python tools for Location Uncertainty in NEMO: read ocean '
            'data, operate on ORCA grids, estimate the statistics and structure of the stochastic noise, reduce records '
            'to a few modes, and diagnose the results.</p><p>Each kit below has a page saying why it exists and what it '
            'holds; each routine has a page with what it does and why, its formula, an example from the tests, where '
            f'the library uses it, and the NEMO lines it follows (searched in NEMO {releases}, newest first).</p><table>'
            + "".join(f'<tr><td><a href="#{k}">{k}</a></td><td>{esc(first_sentence(packages.get(k)))} '
                      f'<span class="toc">{sum(1 for (kk, _) in mods if kk == k)} modules, '
                      f'{sum(len(functions(t)) for (kk, _), t in mods.items() if kk == k)} functions</span></td></tr>' for k in kits)
            + "</table>"
            + (f"<h2>Tests</h2><p>{totals[-1][1]} checks in {totals[-1][0]} test groups of <code>tests/test_all_routines.ipynb</code>, "
               f"as last executed. {len(untested)} routines have no test of their own:</p><p>"
               + ", ".join(f'<a href="#{k}"><code>{esc(k)}</code></a>' for k in untested) + "</p>" if totals else "")
            + "</section>")
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1"><title>py-kits reference</title>'
            f'<link rel="stylesheet" href="{KATEX}/katex.min.css">'
            f'<script defer src="{KATEX}/katex.min.js"></script>'
            f'<script defer src="{KATEX}/contrib/auto-render.min.js" '
            "onload=\"renderMathInElement(document.body, {delimiters: [{left: '$$', right: '$$', display: true}], throwOnError: false})\"></script>"
            f"<style>{STYLE}</style></head><body><nav><a href=\"#home\" style=\"padding-left:0\"><b>py-kits</b></a>{''.join(nav)}</nav>"
            f"<main>{home}{''.join(pages)}</main></body></html>\n")
    dst = os.path.join("docs", "py-kits.html")
    with open(dst, "w") as f:
        f.write(page)
    print(f"{dst}: {count} functions in {len(mods)} modules of {len(kits)} kits, NEMO {releases}")
    for m in nemo.missing:
        print("NEMO line not found:", m, file=sys.stderr)
    sys.exit(1 if nemo.missing else 0)


if __name__ == "__main__":
    main()
