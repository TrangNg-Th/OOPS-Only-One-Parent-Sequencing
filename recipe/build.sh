#!/bin/bash
# Runs inside conda-build's isolated build environment.
# $SRC_DIR  = a copy of the repo 
# $PREFIX   = the environment being built into
set -euo pipefail

# Where the pipeline's files will live inside any environment that installs
# this package. Kept relative to $PREFIX so it works no matter where conda
# puts the environment.
INSTALL_DIR="${PREFIX}/share/oops"

mkdir -p "${INSTALL_DIR}"
cp -r "${SRC_DIR}/main.sh"        "${INSTALL_DIR}/"
cp -r "${SRC_DIR}/src"            "${INSTALL_DIR}/"
cp -r "${SRC_DIR}/data"           "${INSTALL_DIR}/"
cp    "${SRC_DIR}/environment.yml" "${INSTALL_DIR}/"
[ -f "${SRC_DIR}/README.md" ] && cp "${SRC_DIR}/README.md" "${INSTALL_DIR}/"

chmod +x "${INSTALL_DIR}/main.sh"
find "${INSTALL_DIR}/src" -name "*.py" -exec chmod +x {} \;
find "${INSTALL_DIR}/src/slurm_scripts_helper" -name "*.slurm" -exec chmod +x {} \; 2>/dev/null || true

# A thin wrapper so users get a plain `oops` command on PATH instead of
# having to remember where conda put the real script. Resolves its own
# location at runtime, so it survives the environment being moved/cloned.
mkdir -p "${PREFIX}/bin"
cat > "${PREFIX}/bin/oops" <<'WRAPPER'
#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
exec bash "${SCRIPT_DIR}/../share/oops/main.sh" "$@"
WRAPPER
chmod +x "${PREFIX}/bin/oops"
