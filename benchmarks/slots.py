"""Shared slot pool for solver runs (threads and memory), used by run.py and ab.py.

POPSO_SLOTS=<dir>:<threads>:<mem_mb> turns it on: every solver run takes its threads and memory limit
from the pool before it starts (and before its clock starts) and gives them back when it ends, so
several run.py processes (the sides and seeds of an A/B, or several queued jobs) can share one machine
without oversubscribing it. Requests are granted strictly in order of (POPSO_SLOT_PRIO, time of
request): a later job only gets what the earlier jobs' runs leave free, i.e. it fills their tail.

The pool lives in <dir>: a lock file serializing grants, one file per waiting request and per running
holder, each flock'ed by its owner while it lives. A file whose lock can be taken belongs to a process
that died and is removed, so a crash never leaks threads or memory. Within one process only one
thread at a time waits on the pool (the others queue on a thread lock), so a waiting pool costs a few
scans a second in all.
"""

from __future__ import annotations

import fcntl
import itertools
import os
import threading
import time
from contextlib import contextmanager
from pathlib import Path

_seq = itertools.count()
_proc_lock = threading.Lock()  # one waiting thread per process
_POLL = 0.25


def config():
    """(dir, threads, mem_mb) from POPSO_SLOTS, or None."""
    spec = os.environ.get("POPSO_SLOTS")
    if not spec:
        return None
    d, threads, mem = spec.rsplit(":", 2)
    return Path(d), int(threads), int(mem)


_START = time.time()


def _prio() -> float:
    p = os.environ.get("POPSO_SLOT_PRIO")
    return float(p) if p else _START


def _alive(path: Path) -> bool:
    """True if the owner of `path` still holds its lock (false for a missing file)."""
    try:
        fd = os.open(path, os.O_RDWR)
    except FileNotFoundError:
        return False
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return True
    finally:
        os.close(fd)  # also drops the lock if we got it
    return False


def _scan(d: Path):
    """Live waiters [(key, name)] and the threads / memory held by live holders; removes dead files."""
    waiters, used_t, used_m = [], 0, 0
    for p in d.iterdir():
        name = p.name
        if not (name.startswith("w_") or name.startswith("h_")):
            continue
        if not _alive(p):
            try:
                p.unlink()
            except FileNotFoundError:
                pass
            continue
        parts = name.split("_")
        if name.startswith("w_"):  # w_<prio>_<time>_<pid>_<seq>
            waiters.append(((float(parts[1]), float(parts[2]), int(parts[3]), int(parts[4])), name))
        else:  # h_<threads>_<mem>_<pid>_<seq>
            used_t += int(parts[1])
            used_m += int(parts[2])
    waiters.sort()
    return waiters, used_t, used_m


def _own(path: Path) -> int:
    fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
    fcntl.flock(fd, fcntl.LOCK_EX)
    return fd


@contextmanager
def slot(threads: int, mem_mb: int):
    """Holds `threads` threads and `mem_mb` MB of the pool for the duration (no-op without POPSO_SLOTS)."""
    cfg = config()
    if cfg is None:
        yield
        return
    d, tot_t, tot_m = cfg
    d.mkdir(parents=True, exist_ok=True)
    need_t, need_m = min(max(threads, 1), tot_t), min(max(mem_mb, 0), tot_m)
    pid = os.getpid()
    hold_fd = None
    hold_path = None
    with _proc_lock:
        wname = f"w_{_prio():.3f}_{time.time():.6f}_{pid}_{next(_seq)}"
        wpath = d / wname
        # Files are created under the pool lock: a scan must never see one before its owner locks it.
        lfd = os.open(d / "lock", os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(lfd, fcntl.LOCK_EX)
            wfd = _own(wpath)
        finally:
            os.close(lfd)
        try:
            while hold_fd is None:
                lfd = os.open(d / "lock", os.O_RDWR | os.O_CREAT, 0o644)
                try:
                    fcntl.flock(lfd, fcntl.LOCK_EX)
                    waiters, used_t, used_m = _scan(d)
                    head = waiters[0][1] if waiters else None
                    if head == wname and used_t + need_t <= tot_t and used_m + need_m <= tot_m:
                        hold_path = d / f"h_{need_t}_{need_m}_{pid}_{next(_seq)}"
                        hold_fd = _own(hold_path)
                finally:
                    os.close(lfd)
                if hold_fd is None:
                    time.sleep(_POLL)
        finally:
            try:
                wpath.unlink()
            except FileNotFoundError:
                pass
            os.close(wfd)
    try:
        yield
    finally:
        try:
            hold_path.unlink()
        except FileNotFoundError:
            pass
        os.close(hold_fd)
