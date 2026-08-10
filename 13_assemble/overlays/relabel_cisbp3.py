"""Relabel the DBD source string 'CIS-BP 2.0' -> 'CIS-BP 3.0' across the built
per-TF JSONs. CIS-BP was updated to 3.0 with no changes to the DBD calls, so
this is a pure label swap (same-length substring, JSON structure untouched)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tfdir = os.path.join(ROOT, "public/mock/tfs")
OLD, NEW = "CIS-BP 2.0", "CIS-BP 3.0"
changed = 0
for fn in os.listdir(tfdir):
    if not fn.endswith(".json"):
        continue
    p = os.path.join(tfdir, fn)
    txt = open(p, encoding="utf-8").read()
    if OLD in txt:
        open(p, "w", encoding="utf-8").write(txt.replace(OLD, NEW))
        changed += 1
# index.json too, if present
idx = os.path.join(ROOT, "public/mock/index.json")
t = open(idx, encoding="utf-8").read()
if OLD in t:
    open(idx, "w", encoding="utf-8").write(t.replace(OLD, NEW))
    print("index.json: relabeled")
print(f"relabeled CIS-BP 2.0 -> 3.0 in {changed} per-TF JSONs")
