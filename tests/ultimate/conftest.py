# -*- coding: utf-8 -*-
"""
Pytest configuration for the ultimate end-to-end test.

Adds an autouse fixture that fails the ultimate test if any
``logging.CRITICAL`` message is emitted during execution.
"""
from __future__ import annotations

import logging
import pytest


class _CriticalLogCatcher(logging.Handler):
    """Collects every CRITICAL-level log record."""

    def __init__(self):
        super().__init__(level=logging.CRITICAL)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


@pytest.fixture(autouse=True)
def fail_on_critical_log():
    """
    Install a log handler before each test and fail if CRITICAL messages
    were emitted.  This turns silent model-quality problems into explicit
    test failures.
    """
    handler = _CriticalLogCatcher()
    root = logging.getLogger()
    root.addHandler(handler)
    yield
    root.removeHandler(handler)
    if handler.records:
        messages = "\n".join(
            f"  [{r.name}] {r.getMessage()}" for r in handler.records
        )
        pytest.fail(
            f"Test emitted {len(handler.records)} CRITICAL log message(s):\n"
            + messages
        )
