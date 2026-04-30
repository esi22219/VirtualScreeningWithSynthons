from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence
import json
import os
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

# summary to print at end of enumeration if saving to parquet
@dataclass
class OutputSummary:
    out_dir: Optional[str]
    shards: List[str]
    n_records: int
    n_batches: int
    extra: Dict[str, Any]

# sink for saving to parquet files
class ParquetSink:
    def __init__(self, out_dir: str, prefix: str = 'part', compression: str = 'zstd'):
        self.out_dir = str(out_dir)
        self.prefix = prefix
        self.compression = compression
        self._shards: List[str] = []
        self._n_records = 0
        self._n_batches = 0
        Path(self.out_dir).mkdir(parents=True, exist_ok=True)

    def consume(self, batch: List[Dict[str, Any]]) -> None:
        if not batch:
            return
        shard_name = f"{self.prefix}-{self._n_batches:06d}.parquet"
        shard_path = os.path.join(self.out_dir, shard_name)
        table = pa.Table.from_pylist(batch)
        pq.write_table(table, shard_path, compression=self.compression)
        self._shards.append(shard_path)
        self._n_records += len(batch)
        self._n_batches += 1

    def finalize(self, extra: Optional[Dict[str, Any]] = None, write_manifest: bool = True) -> OutputSummary:
        extra = extra or {} # if not exist
        summary = OutputSummary(
            out_dir=self.out_dir,
            shards=list(self._shards),
            n_records=self._n_records,
            n_batches=self._n_batches,
            extra=extra,
        )
        if write_manifest:
            manifest_path = os.path.join(self.out_dir, 'manifest.json')
            with open(manifest_path, 'w') as f:
                json.dump({
                    'out_dir': summary.out_dir,
                    'shards': summary.shards,
                    'n_records': summary.n_records,
                    'n_batches': summary.n_batches,
                    'extra': summary.extra,
                }, f, indent=2)
        return summary
