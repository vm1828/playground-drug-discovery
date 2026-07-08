#!/usr/bin/env bash
# Set up the local environment for the drug-discovery project.
#
#   ./init.sh          # create .venv and install dependencies
#
# PaDEL (via padelpy) needs a Java runtime on your PATH. On Debian/Ubuntu:
#   sudo apt-get install -y default-jre
set -euo pipefail

cd "$(dirname "$0")"

PYTHON="${PYTHON:-python3}"
VENV_DIR="${VENV_DIR:-.venv}"

if [ ! -d "$VENV_DIR" ]; then
  echo "Creating virtual environment in $VENV_DIR ..."
  "$PYTHON" -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

echo "Upgrading pip ..."
python -m pip install --upgrade pip

echo "Installing dependencies ..."
python -m pip install -r requirements.txt

# Register the venv as a Jupyter kernel so VS Code / Jupyter can select it.
python -m ipykernel install --user --name drug-discovery --display-name "Python (drug-discovery)"

if ! command -v java >/dev/null 2>&1; then
  echo
  echo "WARNING: Java was not found on PATH. padelpy needs it to compute PaDEL"
  echo "         fingerprints. Install a JRE, e.g. 'sudo apt-get install default-jre'."
  echo "         (The cached descriptors in data/ let the notebook run without it.)"
fi

echo
echo "Done. Activate with:  source $VENV_DIR/bin/activate"
