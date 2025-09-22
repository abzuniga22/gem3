import os, sqlite3, re, csv, io
from math import ceil
from flask import (
    Flask, render_template, request, abort,
    send_from_directory, url_for, jsonify, Response
)
from markupsafe import Markup, escape

APP_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.getenv("DATA_DIR", APP_DIR)
DB_PATH  = os.path.join(DATA_DIR, "GEM3.db")

print("APP_DIR :", APP_DIR)
print("DATA_DIR:", DATA_DIR)
print("DB_PATH :", DB_PATH)

app = Flask(__name__)

# --- autolink filter ---
# IDs we’ll link ---
# KEGG compound/reaction (allow optional compartment suffix like _c/_m/etc., not part of link)
_KEGG_C = re.compile(r'(?<!\w)(C\d{5})(?:_([a-z]))?(?!\w)', re.I)
_KEGG_R = re.compile(r'(?<!\w)(R\d{5})(?:_([a-z]))?(?!\w)', re.I)

# KEGG pathways: map00010, hsa00010, eco00010
_KEGG_PATH = re.compile(r'(?<!\w)(?:map|hsa|eco)\d{5}(?!\w)', re.I)

# PubChem CID: “CID: 5793”, “PubChem 5793”
_PUBCHEM = re.compile(r'\b(?:CID|PubChem)\s*:?\s*(\d+)\b', re.I)

# ChEBI / HMDB
_CHEBI = re.compile(r'\bCHEBI:(\d+)\b', re.I)
_HMDB  = re.compile(r'\b(HMDB\d{5})\b')

# EC numbers: EC 1.1.1.1 / EC:1.1.1.1
_EC    = re.compile(r'\bEC[:\s]+(\d+\.\d+\.\d+\.\d+)\b', re.I)

# MetaNetX
_MNXM  = re.compile(r'\b(MNXM\d+)\b')  # metabolites
_MNXR  = re.compile(r'\b(MNXR\d+)\b')  # reactions

# ModelSEED
_SEED_C = re.compile(r'\b(cpd\d{5})\b', re.I)
_SEED_R = re.compile(r'\b(rxn\d{5})\b', re.I)

# BiGG metabolites like "glc__D_c", "h2o_c" (lowercase base; single-letter compartment)
_BIGG_M = re.compile(r'\b([a-z0-9_]{2,})_([acgmpexlrns])\b')

# Optional: KEGG gene locus tags like Avin_00080 → avn:Avin_00080
_KEGG_ORG_FOR_LOCUS = {"Avin": "avn"}  # add more mappings when you need them
_LOCUS = re.compile(r'\b([A-Za-z][A-Za-z0-9]{2,10})_([0-9]{3,6})\b')

def autolink_external_ids(value):
    if value is None:
        return ""
    s = escape(str(value))

    # KEGG C/R (keep suffix outside the link if present)
    s = _KEGG_C.sub(lambda m: (
        f'<a target="_blank" rel="noopener" href="https://www.kegg.jp/entry/{m.group(1)}">{m.group(1)}</a>'
        + (f'_{m.group(2)}' if m.group(2) else '')
    ), s)

    s = _KEGG_R.sub(lambda m: (
        f'<a target="_blank" rel="noopener" href="https://www.kegg.jp/entry/{m.group(1)}">{m.group(1)}</a>'
        + (f'_{m.group(2)}' if m.group(2) else '')
    ), s)

    # KEGG pathway
    s = _KEGG_PATH.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://www.kegg.jp/pathway/{m.group(0)}">{m.group(0)}</a>', s)

    # PubChem / ChEBI / HMDB / EC
    s = _PUBCHEM.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://pubchem.ncbi.nlm.nih.gov/compound/{m.group(1)}">{m.group(0)}</a>', s)
    s = _CHEBI.sub(lambda m:
    f'<a target="_blank" rel="noopener" href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{m.group(1)}">CHEBI:{m.group(1)}</a>', s)

    s = _HMDB.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://hmdb.ca/metabolites/{m.group(1)}">{m.group(1)}</a>', s)
    s = _EC.sub(lambda m:
    f'<a target="_blank" rel="noopener" href="https://www.kegg.jp/entry/{m.group(1)}">EC {m.group(1)}</a>', s)

    # MetaNetX & ModelSEED
    s = _MNXM.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://www.metanetx.org/chem_info/{m.group(1)}">{m.group(1)}</a>', s)
    s = _MNXR.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://www.metanetx.org/reac_info/{m.group(1)}">{m.group(1)}</a>', s)
    s = _SEED_C.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://modelseed.org/biochem/compounds/{m.group(1)}">{m.group(1)}</a>', s)
    s = _SEED_R.sub(lambda m:
        f'<a target="_blank" rel="noopener" href="https://modelseed.org/biochem/reactions/{m.group(1)}">{m.group(1)}</a>', s)

    # BiGG metabolite
    s = _BIGG_M.sub(lambda m: (
        f'<a target="_blank" rel="noopener" href="http://bigg.ucsd.edu/universal/metabolites/{m.group(1)}">{m.group(1)}</a>'
        f'_{m.group(2)}'
    ), s)

    # KEGG gene locus (e.g., Avin_00080)
    def _link_locus(m):
        prefix = m.group(1); locus = f"{prefix}_{m.group(2)}"
        org = _KEGG_ORG_FOR_LOCUS.get(prefix)
        url = (f"https://www.kegg.jp/dbget-bin/www_bget?{org}:{locus}"
               if org else f"https://www.kegg.jp/dbget-bin/www_bfind?genes={locus}")
        return f'<a target="_blank" rel="noopener" href="{url}">{locus}</a>'
    s = _LOCUS.sub(_link_locus, s)

    return Markup(s)

app.jinja_env.filters["autolink"] = autolink_external_ids


def get_db():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON;")
    return conn


@app.route("/debug/health")
def debug_health():
    p = DB_PATH
    out = [f"DB_PATH={p}", f"EXISTS={os.path.exists(p)}", f"SIZE={os.path.getsize(p) if os.path.exists(p) else 0}"]
    try:
        with sqlite3.connect(p) as c:
            names = [r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type in ('table','view') LIMIT 10")]
        out.append("TABLES=" + ", ".join(names))
    except Exception as e:
        out.append(f"SQLITE_ERR={e!r}")
    return "<pre>" + "\n".join(out) + "</pre>"


# --- Model → primary paper URL (quick, no-DB solution) ---
MODEL_PAPERS = {
    "CLasA4_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasA4_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasA4_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasA4_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasF17_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasF17_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasF17_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasF17_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasIshi_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasIshi_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasIshi_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasIshi_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasYCPsy_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasYCPsy_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasYCPsy_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasYCPsy_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasgxpsy_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasgxpsy_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasgxpsy_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLasgxpsy_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLaspsy62_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLaspsy62_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLaspsy62_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "CLaspsy62_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "LcBT1_BM7": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "LcBT1_M13": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "LcBT1_M14": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "LcBT1_M15": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7403731/",
    "SA300_a": "https://pubmed.ncbi.nlm.nih.gov/30625152/",
    "SA300_an": "https://pubmed.ncbi.nlm.nih.gov/30625152/",
    "iCZ843_H2": "https://pubmed.ncbi.nlm.nih.gov/27372244/",
    "iCZ843_M2": "https://pubmed.ncbi.nlm.nih.gov/27372244/",
    "iCZ843_PA2": "https://pubmed.ncbi.nlm.nih.gov/27372244/",
    "iDT1278": "https://pubmed.ncbi.nlm.nih.gov/32551229/",
    "iGC535_Benzene": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_Chlorobenzene": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_FructoseHigh": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_FructoseLow": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_HCO3High": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_HCO3Low": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_Phenol": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_Pyruvate": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
    "iGC535_Toluene": "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009828",
}

def paper_url(model_name: str) -> str | None:
    return MODEL_PAPERS.get(model_name)

# expose to Jinja
app.jinja_env.globals["paper_url"] = paper_url

MODELS_DIR = os.path.join(app.root_path, "static", "models")
EXT_ORDER = {".xml": 0, ".mat": 1, ".json": 2}  # desired display order


def get_model_groups():
    try:
        files = [f for f in os.listdir(MODELS_DIR)
                 if os.path.isfile(os.path.join(MODELS_DIR, f))]
    except FileNotFoundError:
        files = []

    groups = {}
    for f in files:
        name, ext = os.path.splitext(f)
        # only include the three types you care about
        if ext.lower() not in EXT_ORDER:
            continue
        groups.setdefault(name, []).append(f)

    # turn into a sorted structure: [(model_name, [(file_name, url), ...]), ...]
    out = []
    for model_name, flist in groups.items():
        flist.sort(key=lambda fn: EXT_ORDER.get(os.path.splitext(fn)[1].lower(), 99))
        out.append((
            model_name,
            [(fn, url_for("download_model", filename=fn)) for fn in flist]
        ))
    out.sort(key=lambda x: x[0].lower())  # sort models alphabetically
    return out


@app.route("/model")
def model_page():
    model_groups = get_model_groups()
    return render_template("model.html", model_groups=model_groups)


@app.route("/download/<path:filename>")
def download_model(filename):
    safe_name = os.path.basename(filename)
    full_path = os.path.join(MODELS_DIR, safe_name)
    if not os.path.isfile(full_path):
        abort(404)
    # send as attachment
    return send_from_directory(
        MODELS_DIR,
        safe_name,
        as_attachment=True,
        download_name=safe_name,
        mimetype="application/octet-stream",
        max_age=31536000
    )


# ------------------ Search (home) ------------------
@app.route("/", methods=["GET"])
def index():
    q = (request.args.get("q") or "").strip()
    scopes = request.args.getlist("scope") or ["model", "reaction", "metabolite", "gene"]

    results = {"model": [], "reaction": [], "metabolite": [], "gene": []}
    if q:
        like = f"%{q}%"
        with get_db() as conn:
            c = conn.cursor()

            if "model" in scopes:
                c.execute(
                    """
                    SELECT model_name, organism, pro_euk, respiration, cell_type
                    FROM model_name
                    WHERE model_name LIKE ? OR IFNULL(organism,'') LIKE ?
                    ORDER BY model_name
                    LIMIT 50
                    """,
                    (like, like),
                )
                results["model"] = c.fetchall()

            if "reaction" in scopes:
                c.execute(
                    """
                    SELECT model_name, reaction_id, abbreviation,
                           IFNULL(reaction,'')   AS equation,
                           IFNULL(pathway,'')    AS pathway,
                           IFNULL(reversible,'') AS reversible
                    FROM reaction
                    WHERE reaction_id LIKE ?
                       OR IFNULL(reaction,'') LIKE ?
                       OR IFNULL(pathway,'')  LIKE ?
                    ORDER BY model_name, reaction_id
                    LIMIT 50
                    """,
                    (like, like, like),
                )
                results["reaction"] = c.fetchall()

            if "metabolite" in scopes:
                c.execute(
                    """
                    SELECT model_name, metabolite_id, name,
                           IFNULL(formula,'')     AS formula,
                           IFNULL(compartment,'') AS compartment
                    FROM metabolite
                    WHERE metabolite_id LIKE ?
                       OR name LIKE ?
                       OR IFNULL(formula,'') LIKE ?
                    ORDER BY model_name, metabolite_id
                    LIMIT 50
                    """,
                    (like, like, like),
                )
                results["metabolite"] = c.fetchall()

            if "gene" in scopes:
                c.execute(
                    """
                    SELECT model_name, gene_id, IFNULL(protein,'') AS protein
                    FROM gene
                    WHERE gene_id LIKE ?
                       OR IFNULL(protein,'') LIKE ?
                    ORDER BY model_name, gene_id
                    LIMIT 50
                    """,
                    (like, like),
                )
                results["gene"] = c.fetchall()

    return render_template("index.html", q=q, scopes=scopes, results=results)


# ------------------ Model detail with pagination ------------------
@app.route("/model/<model>")
def model_detail(model):
    PAGE_SIZE = 25

    def get_page(arg):
        try:
            p = int(request.args.get(arg, 1))
            return p if p > 0 else 1
        except (TypeError, ValueError):
            return 1

    rpage = get_page("rpage")
    mpage = get_page("mpage")
    gpage = get_page("gpage")

    with get_db() as conn:
        c = conn.cursor()

        # Model info
        m = c.execute(
            """
            SELECT model_name AS model,
                   COALESCE(organism,'')     AS organism,
                   COALESCE(pro_euk,'')      AS pro_euk,
                   COALESCE(respiration,'')  AS respiration,
                   COALESCE(cell_type,'')    AS cell_type
            FROM model_name
            WHERE model_name = ?
            """,
            (model,),
        ).fetchone()
        if not m:
            abort(404)

        # Counts
        n_rxn = c.execute("SELECT COUNT(*) AS n FROM reaction   WHERE model_name=?", (model,)).fetchone()["n"]
        n_met = c.execute("SELECT COUNT(*) AS n FROM metabolite WHERE model_name=?", (model,)).fetchone()["n"]
        n_gene = c.execute("SELECT COUNT(*) AS n FROM gene      WHERE model_name=?", (model,)).fetchone()["n"]

        # Total pages
        rpages = max(ceil(n_rxn / PAGE_SIZE), 1)
        mpages = max(ceil(n_met / PAGE_SIZE), 1)
        gpages = max(ceil(n_gene / PAGE_SIZE), 1)

        # Clamp
        rpage = min(rpage, rpages)
        mpage = min(mpage, mpages)
        gpage = min(gpage, gpages)

        # Offsets
        r_off = (rpage - 1) * PAGE_SIZE
        m_off = (mpage - 1) * PAGE_SIZE
        g_off = (gpage - 1) * PAGE_SIZE

        # Slices
        reactions = c.execute(
            """
            SELECT reaction_id,
                   IFNULL(reaction,'') AS equation,
                   IFNULL(pathway,'')  AS pathway
            FROM reaction
            WHERE model_name = ?
            ORDER BY reaction_id
            LIMIT ? OFFSET ?
            """,
            (model, PAGE_SIZE, r_off),
        ).fetchall()

        metabolites = c.execute(
            """
            SELECT metabolite_id, name,
                   IFNULL(formula,'')     AS formula,
                   IFNULL(compartment,'') AS compartment
            FROM metabolite
            WHERE model_name = ?
            ORDER BY metabolite_id
            LIMIT ? OFFSET ?
            """,
            (model, PAGE_SIZE, m_off),
        ).fetchall()

        genes = c.execute(
            """
            SELECT gene_id, IFNULL(protein,'') AS protein
            FROM gene
            WHERE model_name = ?
            ORDER BY gene_id
            LIMIT ? OFFSET ?
            """,
            (model, PAGE_SIZE, g_off),
        ).fetchall()

    return render_template(
        "model.html",
        m=m,
        n_rxn=n_rxn, n_met=n_met, n_gene=n_gene,
        reactions=reactions, metabolites=metabolites, genes=genes,
        rpage=rpage, rpages=rpages,
        mpage=mpage, mpages=mpages,
        gpage=gpage, gpages=gpages,
    )


# ------------------ NEW: lightweight JSON API ------------------
# GET /api/search?q=...&scope=model&scope=reaction&limit=50&offset=0
@app.route("/api/search")
def api_search():
    q = (request.args.get("q") or "").strip()
    scopes = request.args.getlist("scope") or ["model", "reaction", "metabolite", "gene"]
    limit = max(1, min(int(request.args.get("limit", 50)), 200))
    offset = max(0, int(request.args.get("offset", 0)))

    data = {"model": [], "reaction": [], "metabolite": [], "gene": []}
    total = {"model": 0, "reaction": 0, "metabolite": 0, "gene": 0}

    if not q:
        return jsonify({"q": q, "scopes": scopes, "total": total, "results": data})

    like = f"%{q}%"
    with get_db() as conn:
        c = conn.cursor()

        if "model" in scopes:
            total["model"] = c.execute(
                "SELECT COUNT(*) FROM model_name WHERE model_name LIKE ? OR IFNULL(organism,'') LIKE ?",
                (like, like),
            ).fetchone()[0]
            data["model"] = [dict(r) for r in c.execute(
                """
                SELECT model_name, organism, pro_euk, respiration, cell_type
                FROM model_name
                WHERE model_name LIKE ? OR IFNULL(organism,'') LIKE ?
                ORDER BY model_name
                LIMIT ? OFFSET ?
                """,
                (like, like, limit, offset),
            ).fetchall()]

        if "reaction" in scopes:
            total["reaction"] = c.execute(
                """
                SELECT COUNT(*) FROM reaction
                WHERE reaction_id LIKE ? OR IFNULL(reaction,'') LIKE ? OR IFNULL(pathway,'') LIKE ?
                """,
                (like, like, like),
            ).fetchone()[0]
            data["reaction"] = [dict(r) for r in c.execute(
                """
                SELECT model_name, reaction_id, abbreviation,
                       IFNULL(reaction,'')   AS equation,
                       IFNULL(pathway,'')    AS pathway,
                       IFNULL(reversible,'') AS reversible
                FROM reaction
                WHERE reaction_id LIKE ? OR IFNULL(reaction,'') LIKE ? OR IFNULL(pathway,'') LIKE ?
                ORDER BY model_name, reaction_id
                LIMIT ? OFFSET ?
                """,
                (like, like, like, limit, offset),
            ).fetchall()]

        if "metabolite" in scopes:
            total["metabolite"] = c.execute(
                """
                SELECT COUNT(*) FROM metabolite
                WHERE metabolite_id LIKE ? OR name LIKE ? OR IFNULL(formula,'') LIKE ?
                """,
                (like, like, like),
            ).fetchone()[0]
            data["metabolite"] = [dict(r) for r in c.execute(
                """
                SELECT model_name, metabolite_id, name,
                       IFNULL(formula,'')     AS formula,
                       IFNULL(compartment,'') AS compartment
                FROM metabolite
                WHERE metabolite_id LIKE ? OR name LIKE ? OR IFNULL(formula,'') LIKE ?
                ORDER BY model_name, metabolite_id
                LIMIT ? OFFSET ?
                """,
                (like, like, like, limit, offset),
            ).fetchall()]

        if "gene" in scopes:
            total["gene"] = c.execute(
                "SELECT COUNT(*) FROM gene WHERE gene_id LIKE ? OR IFNULL(protein,'') LIKE ?",
                (like, like),
            ).fetchone()[0]
            data["gene"] = [dict(r) for r in c.execute(
                """
                SELECT model_name, gene_id, IFNULL(protein,'') AS protein
                FROM gene
                WHERE gene_id LIKE ? OR IFNULL(protein,'') LIKE ?
                ORDER BY model_name, gene_id
                LIMIT ? OFFSET ?
                """,
                (like, like, limit, offset),
            ).fetchall()]

    return jsonify({"q": q, "scopes": scopes, "total": total, "results": data})


# GET /api/model/<model>.json?rpage=1&mpage=1&gpage=1
@app.route("/api/model/<model>.json")
def api_model(model):
    PAGE_SIZE = int(request.args.get("page_size", 25))

    def get_page(arg):
        try:
            p = int(request.args.get(arg, 1))
            return p if p > 0 else 1
        except (TypeError, ValueError):
            return 1

    rpage = get_page("rpage")
    mpage = get_page("mpage")
    gpage = get_page("gpage")

    with get_db() as conn:
        c = conn.cursor()

        m = c.execute(
            """
            SELECT model_name AS model,
                   COALESCE(organism,'')     AS organism,
                   COALESCE(pro_euk,'')      AS pro_euk,
                   COALESCE(respiration,'')  AS respiration,
                   COALESCE(cell_type,'')    AS cell_type
            FROM model_name
            WHERE model_name = ?
            """,
            (model,),
        ).fetchone()
        if not m:
            abort(404)

        n_rxn = c.execute("SELECT COUNT(*) AS n FROM reaction   WHERE model_name=?", (model,)).fetchone()["n"]
        n_met = c.execute("SELECT COUNT(*) AS n FROM metabolite WHERE model_name=?", (model,)).fetchone()["n"]
        n_gene = c.execute("SELECT COUNT(*) AS n FROM gene      WHERE model_name=?", (model,)).fetchone()["n"]

        def clamp(page, n):
            pages = max(ceil(n / PAGE_SIZE), 1)
            return min(page, pages), pages

        rpage, rpages = clamp(rpage, n_rxn)
        mpage, mpages = clamp(mpage, n_met)
        gpage, gpages = clamp(gpage, n_gene)

        r_off = (rpage - 1) * PAGE_SIZE
        m_off = (mpage - 1) * PAGE_SIZE
        g_off = (gpage - 1) * PAGE_SIZE

        reactions = [dict(r) for r in c.execute(
            """
            SELECT reaction_id,
                   IFNULL(reaction,'') AS equation,
                   IFNULL(pathway,'')  AS pathway
            FROM reaction
            WHERE model_name = ?
            ORDER BY reaction_id
            LIMIT ? OFFSET ?
            """,
            (model, PAGE_SIZE, r_off),
        ).fetchall()]

        metabolites = [dict(r) for r in c.execute(
            """
            SELECT metabolite_id, name,
                   IFNULL(formula,'')     AS formula,
                   IFNULL(compartment,'') AS compartment
            FROM metabolite
            WHERE model_name = ?
            ORDER BY metabolite_id
            LIMIT ? OFFSET ?
            """,
            (model, PAGE_SIZE, m_off),
        ).fetchall()]

        genes = [dict(r) for r in c.execute(
            """
            SELECT gene_id, IFNULL(protein,'') AS protein
            FROM gene
            WHERE model_name = ?
            ORDER BY gene_id
            LIMIT ? OFFSET ?
            """,
            (model, PAGE_SIZE, g_off),
        ).fetchall()]

    return jsonify({
        "model": dict(m),
        "counts": {"reaction": n_rxn, "metabolite": n_met, "gene": n_gene},
        "page_size": PAGE_SIZE,
        "pages": {"rpage": rpage, "rpages": rpages, "mpage": mpage, "mpages": mpages, "gpage": gpage, "gpages": gpages},
        "reactions": reactions,
        "metabolites": metabolites,
        "genes": genes,
    })


# ------------------ NEW: CSV export ------------------
# /export/model.csv?q=...
# /export/reaction.csv?q=...
# /export/metabolite.csv?q=...
# /export/gene.csv?q=...
ALLOWED_SCOPES = {"model", "reaction", "metabolite", "gene"}

@app.route("/export/<scope>.csv")
def export_scope_csv(scope):
    if scope not in ALLOWED_SCOPES:
        abort(404)

    q = (request.args.get("q") or "").strip()
    like = f"%{q}%" if q else None

    def stream_rows():
        buf = io.StringIO()
        w = csv.writer(buf)
        with get_db() as conn:
            c = conn.cursor()
            if scope == "model":
                w.writerow(["model_name", "organism", "pro_euk", "respiration", "cell_type"]) ; yield buf.getvalue(); buf.seek(0); buf.truncate(0)
                sql = "SELECT model_name, organism, pro_euk, respiration, cell_type FROM model_name"
                params = ()
                if like:
                    sql += " WHERE model_name LIKE ? OR IFNULL(organism,'') LIKE ?"; params = (like, like)
                sql += " ORDER BY model_name"
                for row in c.execute(sql, params):
                    w.writerow([row["model_name"], row["organism"], row["pro_euk"], row["respiration"], row["cell_type"]])
                    yield buf.getvalue(); buf.seek(0); buf.truncate(0)

            elif scope == "reaction":
                w.writerow(["model_name", "reaction_id", "abbreviation", "equation", "pathway", "reversible"]) ; yield buf.getvalue(); buf.seek(0); buf.truncate(0)
                sql = (
                    "SELECT model_name, reaction_id, abbreviation, IFNULL(reaction,'') AS equation, "
                    "IFNULL(pathway,'') AS pathway, IFNULL(reversible,'') AS reversible FROM reaction"
                )
                params = ()
                if like:
                    sql += " WHERE reaction_id LIKE ? OR IFNULL(reaction,'') LIKE ? OR IFNULL(pathway,'') LIKE ?"; params = (like, like, like)
                sql += " ORDER BY model_name, reaction_id"
                for row in c.execute(sql, params):
                    w.writerow([row["model_name"], row["reaction_id"], row["abbreviation"], row["equation"], row["pathway"], row["reversible"]])
                    yield buf.getvalue(); buf.seek(0); buf.truncate(0)

            elif scope == "metabolite":
                w.writerow(["model_name", "metabolite_id", "name", "formula", "compartment"]) ; yield buf.getvalue(); buf.seek(0); buf.truncate(0)
                sql = (
                    "SELECT model_name, metabolite_id, name, IFNULL(formula,'') AS formula, IFNULL(compartment,'') AS compartment FROM metabolite"
                )
                params = ()
                if like:
                    sql += " WHERE metabolite_id LIKE ? OR name LIKE ? OR IFNULL(formula,'') LIKE ?"; params = (like, like, like)
                sql += " ORDER BY model_name, metabolite_id"
                for row in c.execute(sql, params):
                    w.writerow([row["model_name"], row["metabolite_id"], row["name"], row["formula"], row["compartment"]])
                    yield buf.getvalue(); buf.seek(0); buf.truncate(0)

            elif scope == "gene":
                w.writerow(["model_name", "gene_id", "protein"]) ; yield buf.getvalue(); buf.seek(0); buf.truncate(0)
                sql = "SELECT model_name, gene_id, IFNULL(protein,'') AS protein FROM gene"
                params = ()
                if like:
                    sql += " WHERE gene_id LIKE ? OR IFNULL(protein,'') LIKE ?"; params = (like, like)
                sql += " ORDER BY model_name, gene_id"
                for row in c.execute(sql, params):
                    w.writerow([row["model_name"], row["gene_id"], row["protein"]])
                    yield buf.getvalue(); buf.seek(0); buf.truncate(0)

    filename = f"{scope}_export.csv"
    headers = {
        "Content-Type": "text/csv; charset=utf-8",
        "Content-Disposition": f"attachment; filename={filename}",
        "Cache-Control": "no-store",
    }
    return Response(stream_rows(), headers=headers)


# ------------------ Local dev entrypoint ------------------
# On Render, prefer: Start Command = `gunicorn app:app` and do NOT call app.run()
if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(os.environ.get("PORT", 5000)), debug=True)
