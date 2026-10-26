"""Solve environment.yml for linux-64 and write the explicit lock that the
Dockerfile installs from. Runs on any machine with conda (no Docker needed).
Rerun after changing environment.yml, then rebuild the image.

    python containers/update_lock.py
"""
import json
import os
import subprocess
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
specs, deps = [], False
for line in open(os.path.join(ROOT, "environment.yml")):
    if line.startswith("dependencies:"):
        deps = True
    elif deps and line.startswith("  - "):
        specs.append(line[4:].strip())
# an empty package cache, so every package is listed under FETCH with its url and md5
cache = tempfile.mkdtemp()
env = dict(os.environ, CONDA_SUBDIR="linux-64", CONDA_OVERRIDE_GLIBC="2.28", CONDA_PKGS_DIRS=cache)
out = subprocess.run(["conda", "create", "-n", "lock", "--dry-run", "--json", "--override-channels",
                      "-c", "conda-forge", "-c", "bioconda"] + specs, env=env, capture_output=True, text=True)
actions = json.loads(out.stdout)["actions"]
fetch = {p["name"]: p for p in actions["FETCH"]}
assert len(fetch) == len(actions["LINK"]), "package cache not empty"
links = [fetch[p["name"]] for p in actions["LINK"]]  # install order
with open(os.path.join(ROOT, "containers/linux-64.lock"), "w") as fh:
    fh.write("# linux-64 explicit lock for environment.yml, written by containers/update_lock.py\n@EXPLICIT\n")
    for p in links:
        fh.write(f"{p['url']}#{p['md5']}\n")
print(len(links), "packages")
