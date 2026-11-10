#!/bin/sh
# Capture golden outputs: build the legacy pipeline from a git ref, on frozen inputs, in the ODK
# container, and save what it produces.
#
#   golden/capture.sh [OUT_DIR]          (default: golden/out)
#
# Environment:
#   GOLDEN_REF   git ref holding the legacy pipeline        (default: main)
#   OMIM_DATA    directory of OMIM's licensed .txt downloads (default: omim-data)
#   ODK_TAG      obolibrary/odkfull tag                      (default: v1.5.3, as in CI)
#
# Inputs come from three places:
#   - OMIM's downloads, from $OMIM_DATA. Licensed, never committed.
#   - Mondo's OMIM mappings, downloaded once into golden/inputs/ and reused, so that later captures
#     and the rewrite's verification run see the same file.
#   - Everything else (curated TSVs, HGNC, the PubMed/mappings cache) from $GOLDEN_REF itself.
set -eu

cd "$(git rev-parse --show-toplevel)"

out=${1:-golden/out}
ref=${GOLDEN_REF:-main}
omim_data=${OMIM_DATA:-omim-data}
image=obolibrary/odkfull:${ODK_TAG:-v1.5.3}
inputs=golden/inputs
work=golden/work

omim_files="mimTitles.txt mim2gene.txt genemap2.txt morbidmap.txt phenotypicSeries.txt"
mondo_files="mondo_exactmatch_omim.sssom.tsv mondo_exactmatch_omimps.sssom.tsv"
# The released artefacts. All are byte-reproducible across builds; the intermediates omim.ttl (random
# blank-node labels) and mondo-omim-genes.tsv (unordered query results) are not, so are not kept.
outputs="omim.owl mondo-omim-genes.robot.tsv disease-gene-relationships-qc.tsv \
mondo-omim-susceptibility-subset.robot.tsv review.tsv"

for f in $omim_files; do
    [ -s "$omim_data/$f" ] || { echo "missing OMIM download: $omim_data/$f" >&2; exit 1; }
done

mkdir -p "$inputs"
for f in $mondo_files; do
    [ -s "$inputs/$f" ] || curl -fsSL "http://purl.obolibrary.org/obo/mondo/mappings/$f" -o "$inputs/$f"
done

rm -rf "$work"
mkdir -p "$work"
git archive "$ref" | tar -x -C "$work"
for f in $omim_files; do cp "$omim_data/$f" "$work/data/"; done
for f in $mondo_files; do cp "$inputs/$f" "$work/"; done
cp golden/run_legacy.py "$work/"

# omim.ttl is built by run_legacy.py rather than by its make target, which re-downloads every input.
# With it in place, make builds only what is derived from it.
docker run --rm -v "$PWD/$work:/work" -w /work -e ROBOT_JAVA_ARGS=-Xmx8G "$image" sh -ec '
    python3 run_legacy.py
    make omim.owl mondo-omim-genes.robot.tsv disease-gene-relationships-qc.tsv
    {
        python3 --version
        robot --version
        pip freeze 2>/dev/null | grep -iE "^(rdflib|pandas|sssom|pyyaml|python-dotenv)=="
    } > versions.txt
'

rm -rf "$out"
mkdir -p "$out"
for f in $outputs; do cp "$work/$f" "$out/"; done

{
    echo "ref: $(git rev-parse "$ref")"
    echo "captured: $(date +%Y-%m-%d)"
    echo "image: $image"
    sed 's/^/tool: /' "$work/versions.txt"
    for f in $omim_files; do echo "input: $(shasum -a 256 "$omim_data/$f" | cut -d' ' -f1)  $f"; done
    for f in $mondo_files; do echo "input: $(shasum -a 256 "$inputs/$f" | cut -d' ' -f1)  $f"; done
    for f in $outputs; do echo "output: $(shasum -a 256 "$out/$f" | cut -d' ' -f1)  $f"; done
} > "$out/MANIFEST"

rm -rf "$work"
echo "golden outputs written to $out"
