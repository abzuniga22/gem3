
import os
import sqlite3
from math import ceil
from flask import Flask, render_template, request, abort




APP_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.getenv("DATA_DIR", APP_DIR)  # we'll set DATA_DIR on Render
DB_PATH  = os.path.join(DATA_DIR, "GEM3.db")
app = Flask(__name__)

def get_db():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON;")
    return conn

@app.route("/debug/health")
def debug_health():
    import os, sqlite3
    p = DB_PATH
    out = [f"DB_PATH={p}",
           f"EXISTS={os.path.exists(p)}",
           f"SIZE={os.path.getsize(p) if os.path.exists(p) else 0}"]
    try:
        with sqlite3.connect(p) as c:
            names = [r[0] for r in c.execute(
                "SELECT name FROM sqlite_master WHERE type in ('table','view') LIMIT 10")]
        out.append("TABLES=" + ", ".join(names))
    except Exception as e:
        out.append(f"SQLITE_ERR={e!r}")
    return "<pre>" + "\n".join(out) + "</pre>"

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
                c.execute("""
                    SELECT model_name, organism, pro_euk, respiration, cell_type
                    FROM model_name
                    WHERE model_name LIKE ? OR IFNULL(organism,'') LIKE ?
                    ORDER BY model_name
                    LIMIT 50
                """, (like, like))
                results["model"] = c.fetchall()

            if "reaction" in scopes:
                c.execute("""
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
                """, (like, like, like))
                results["reaction"] = c.fetchall()

            if "metabolite" in scopes:
                c.execute("""
                    SELECT model_name, metabolite_id, name,
                           IFNULL(formula,'')     AS formula,
                           IFNULL(compartment,'') AS compartment
                    FROM metabolite
                    WHERE metabolite_id LIKE ?
                       OR name LIKE ?
                       OR IFNULL(formula,'') LIKE ?
                    ORDER BY model_name, metabolite_id
                    LIMIT 50
                """, (like, like, like))
                results["metabolite"] = c.fetchall()

            if "gene" in scopes:
                c.execute("""
                    SELECT model_name, gene_id, IFNULL(protein,'') AS protein
                    FROM gene
                    WHERE gene_id LIKE ?
                       OR IFNULL(protein,'') LIKE ?
                    ORDER BY model_name, gene_id
                    LIMIT 50
                """, (like, like))
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
        m = c.execute("""
            SELECT model_name AS model,
                   COALESCE(organism,'')     AS organism,
                   COALESCE(pro_euk,'')      AS pro_euk,
                   COALESCE(respiration,'')  AS respiration,
                   COALESCE(cell_type,'')    AS cell_type
            FROM model_name
            WHERE model_name = ?
        """, (model,)).fetchone()
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
        reactions = c.execute("""
            SELECT reaction_id,
                   IFNULL(reaction,'') AS equation,
                   IFNULL(pathway,'')  AS pathway
            FROM reaction
            WHERE model_name = ?
            ORDER BY reaction_id
            LIMIT ? OFFSET ?
        """, (model, PAGE_SIZE, r_off)).fetchall()

        metabolites = c.execute("""
            SELECT metabolite_id, name,
                   IFNULL(formula,'')     AS formula,
                   IFNULL(compartment,'') AS compartment
            FROM metabolite
            WHERE model_name = ?
            ORDER BY metabolite_id
            LIMIT ? OFFSET ?
        """, (model, PAGE_SIZE, m_off)).fetchall()

        genes = c.execute("""
            SELECT gene_id, IFNULL(protein,'') AS protein
            FROM gene
            WHERE model_name = ?
            ORDER BY gene_id
            LIMIT ? OFFSET ?
        """, (model, PAGE_SIZE, g_off)).fetchall()

    return render_template(
        "model.html",
        m=m,
        n_rxn=n_rxn, n_met=n_met, n_gene=n_gene,
        reactions=reactions, metabolites=metabolites, genes=genes,
        rpage=rpage, rpages=rpages,
        mpage=mpage, mpages=mpages,
        gpage=gpage, gpages=gpages
    )

if __name__ == "__main__":
    app.run(debug=True)
