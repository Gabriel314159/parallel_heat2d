#!/usr/bin/env python3
"""Read-only probe executed once on each allocated compute node."""
import json
import os
import socket
from pathlib import Path

affinity = sorted(os.sched_getaffinity(0))
rows = []
for cpu in affinity:
    topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    rows.append({"cpu": cpu,
                 "socket": int((topology / "physical_package_id").read_text()),
                 "core": int((topology / "core_id").read_text())})
expected = set(range(72))
by_cpu = {r["cpu"]: r for r in rows}
ok = expected.issubset(by_cpu)
if ok:
    even = {by_cpu[c]["socket"] for c in range(0, 72, 2)}
    odd = {by_cpu[c]["socket"] for c in range(1, 72, 2)}
    ok = (len(even) == len(odd) == 1 and even != odd and
          len({(by_cpu[c]["socket"], by_cpu[c]["core"]) for c in expected}) == 72)
print(json.dumps({"host": socket.gethostname(), "affinity": affinity,
                  "topology": rows, "vanda_layout_ok": ok}), flush=True)
