import importlib.util
import sys
from pathlib import Path

patched = Path(__file__).with_name("metagene.p95fix.py")

spec = importlib.util.spec_from_file_location(
    "plastid.bin.metagene_p95fix",
    str(patched)
)

mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

ret = mod.main(sys.argv[1:])

if ret is None:
    ret = 0

raise SystemExit(ret)
