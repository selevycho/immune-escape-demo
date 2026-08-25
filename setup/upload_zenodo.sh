#!/bin/bash
#
# Upload the demo package to Zenodo from the command line.
#
# The package lives on a cluster with no desktop, so the web form is not
# an option. Zenodo's REST API takes the same files over HTTPS.
#
# The token is read from the environment and never written anywhere: it
# grants write access to every record the account owns, so it does not
# belong in a file, a log, or a shell history.
#
# Files are uploaded one at a time and each is verified against the
# checksum Zenodo computes on receipt, because a truncated 300 MB BAM
# would otherwise be discovered only when someone tried to run the demo.
#
# Publishing is deliberately not automated. Once a record is published its
# files are immutable, and that decision should be made by a person who
# has seen the upload succeed.
#
# Usage:
#   export ZENODO_TOKEN=...
#   ./upload_zenodo.sh
#   ./upload_zenodo.sh --dry-run
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

PKG="${PKG:-${DATA_DIR}/zenodo}"
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

if [ -z "${ZENODO_TOKEN:-}" ]; then
    fail "ZENODO_TOKEN is not set"
    echo ""
    echo "  create one at"
    echo "    https://zenodo.org/account/settings/applications/tokens/new/"
    echo "  with scopes deposit:write and deposit:actions, then"
    echo ""
    echo "    read -s -p 'token: ' ZENODO_TOKEN; export ZENODO_TOKEN"
    exit 1
fi

if [ ! -d "${PKG}" ]; then
    fail "no package at ${PKG}; run make_zenodo_package.sh first"
    exit 1
fi

API="https://zenodo.org/api"

say "the record"
resp="$(curl -s "${API}/deposit/depositions/${ZENODO_RECORD}?access_token=${ZENODO_TOKEN}")"

if echo "${resp}" | grep -q '"status": *4'; then
    fail "cannot read record ${ZENODO_RECORD}"
    echo "${resp}" | head -c 400
    echo ""
    echo ""
    echo "  the token may lack deposit:write, or the record id in"
    echo "  config.sh may be wrong. ZENODO_RECORD is currently"
    echo "  ${ZENODO_RECORD}, taken from ${ZENODO_DOI}"
    exit 1
fi

PY=""
for cand in python3 "$(conda info --base 2>/dev/null)/envs/${ENV_GATK}/bin/python"; do
    "${cand}" -c "import json" >/dev/null 2>&1 && PY="${cand}" && break
done

# a heredoc consumes stdin, so the response cannot also be piped in;
# grep and sed read the two fields needed without a JSON parser
show() {
    local key="$1" label="$2"
    local v
    v="$(echo "${resp}" | grep -o "\"${key}\": *\"[^\"]*\"" | head -1 \
         | sed 's/.*: *"//; s/"$//')"
    printf '  %-11s%s\n' "${label}" "${v:-(none)}"
}

show title title
show state state
n_files="$(echo "${resp}" | grep -o '"key": *"[^"]*"' | wc -l)"
printf '  %-11s%s already there\n' "files" "${n_files}"

BUCKET="$(echo "${resp}" | grep -o '"bucket": *"[^"]*"' | head -1 \
          | sed 's/.*: *"//; s/"$//')"
printf '  %-11s%s\n' "bucket" "$([ -n "${BUCKET}" ] && echo yes || echo no)"

if [ -z "${BUCKET}" ]; then
    fail "no upload bucket; the record may already be published"
    exit 1
fi

say "files to send"
total=0
for f in "${PKG}"/*; do
    [ -f "${f}" ] || continue
    sz=$(stat -Lc%s "${f}")
    total=$((total + sz))
    printf '  %-28s %8s\n' "$(basename "${f}")" \
        "$(du -Lh "${f}" | cut -f1)"
done
note ""
note "total $(echo "${total}" | awk '{printf "%.1f GB", $1/1073741824}')"

if [ "${DRY}" = "1" ]; then
    echo ""
    echo "  dry run — nothing sent"
    exit 0
fi

say "uploading"
SENT=0; FAILED=0
for f in "${PKG}"/*; do
    [ -f "${f}" ] || continue
    name="$(basename "${f}")"
    printf '  %-28s ' "${name}"

    # the bucket API takes a plain PUT of the file body, which streams
    # rather than buffering the whole thing in memory as a multipart form
    # would
    out="$(curl -s --retry 3 --retry-delay 5 \
             -X PUT \
             -H "Authorization: Bearer ${ZENODO_TOKEN}" \
             --upload-file "${f}" \
             "${BUCKET}/${name}")"

    got="$(echo "${out}" | grep -o '"checksum": *"[^"]*"' | head -1 \
           | sed 's/.*md5://; s/"$//')"

    want="$(md5sum "${f}" | cut -d" " -f1)"

    if [ -n "${got}" ] && [ "${got}" = "${want}" ]; then
        echo "ok"
        SENT=$((SENT+1))
    elif [ -n "${got}" ]; then
        echo "CHECKSUM MISMATCH"
        FAILED=$((FAILED+1))
    else
        echo "failed"
        echo "${out}" | head -c 200
        echo ""
        FAILED=$((FAILED+1))
    fi
done

say "result"
note "${SENT} uploaded, ${FAILED} failed"

if [ "${FAILED}" -gt 0 ]; then
    echo ""
    echo "  rerun to retry; files already there are overwritten cleanly"
    exit 1
fi

echo ""
echo "  All files are on the record and verified."
echo ""
echo "  It is still a draft. Nothing is public and the files can still"
echo "  be replaced. Run the demo end to end first; publish only once"
echo "  you know these are the files it needs."
echo ""
echo "  To publish:"
echo "    https://zenodo.org/uploads/${ZENODO_RECORD}"
echo "  or"
echo "    curl -X POST \\"
echo "      \"${API}/deposit/depositions/${ZENODO_RECORD}/actions/publish?access_token=\$ZENODO_TOKEN\""
