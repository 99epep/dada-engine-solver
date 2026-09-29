"""Durable files with an append-only evaluation journal and recoverable snapshots."""
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import fcntl
import json
import os
from dada_solver.campaign.candidate import canonical_json, content_hash


def atomic_json(path, value):
    atomic_text(path, canonical_json(value)+'\n')


def atomic_text(path, text):
    """Durably publish a text snapshot using the journal's atomic-write contract."""
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    with temporary.open('w') as stream:
        stream.write(text); stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)
    # Persist the rename as well as file contents on local POSIX filesystems.
    fd = os.open(path.parent, os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


class CampaignHistory:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.path = self.directory/'history.jsonl'

    @contextmanager
    def locked(self):
        with (self.directory/'.writer.lock').open('a') as stream:
            try: fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError('Another writer is already running this campaign.') from error
            try: yield
            finally: fcntl.flock(stream, fcntl.LOCK_UN)

    def load(self):
        records = []
        if self.path.exists():
            data = self.path.read_bytes(); offset = 0
            for line in data.splitlines(keepends=True):
                try:
                    record = json.loads(line)
                    if not line.endswith(b'\n'): raise ValueError('Incomplete final journal line.')
                    records.append(record); offset += len(line)
                except (ValueError, UnicodeDecodeError):
                    if offset+len(line) != len(data):
                        raise ValueError('Corrupt non-final history record; manual recovery is required.')
                    # Sole append-only exception: preserve then remove a torn final write.
                    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
                    (self.directory/f'history_torn_tail_{stamp}.bin').write_bytes(data[offset:])
                    with self.path.open('r+b') as stream:
                        stream.truncate(offset); stream.flush(); os.fsync(stream.fileno())
        known = {r['evaluation_number']:r for r in records}
        orphans = []
        for record in recovery_records(self.directory):
            verify_record(record)
            number = record['evaluation_number']
            if number in known:
                if known[number] != record:
                    raise ValueError('Recovery conflicts with a journal evaluation.')
            else:
                orphans.append(record); known[number] = record
        # Validate everything before publishing recovered records.
        numbers = sorted(known)
        if numbers != list(range(len(numbers))):
            raise ValueError('History evaluation numbers are not contiguous.')
        for record in sorted(orphans, key=lambda r:r['evaluation_number']):
            self._append(record); records.append(record)
        records.sort(key=lambda r:r['evaluation_number'])
        for record in records: verify_record(record)
        if [r['evaluation_number'] for r in records] != list(range(len(records))):
            raise ValueError('History evaluation numbers are not contiguous.')
        return records

    def _append(self, record):
        with self.path.open('a') as stream:
            stream.write(canonical_json(record)+'\n'); stream.flush(); os.fsync(stream.fileno())

    def save(self, record):
        # One recovery slot, including cache hits: durable completion precedes
        # the append. The runner acknowledges only after state.json is durable.
        atomic_json(self.directory/'recovery.json', record)
        self._append(record)

    def clear_recovery(self):
        path = self.directory/'recovery.json'
        if path.exists():
            path.unlink()
            fd = os.open(self.directory, os.O_RDONLY)
            try: os.fsync(fd)
            finally: os.close(fd)

    def state(self):
        path = self.directory/'state.json'
        return json.loads(path.read_text()) if path.exists() else {}

    def save_state(self, state):
        atomic_json(self.directory/'state.json', state)


def verify_record(record):
    payload = {k:record[k] for k in ('schema_version','definition_id','normalized',
        'physical','families','numerical_settings')}
    if content_hash(payload) != record['candidate_id']:
        raise ValueError('Persisted candidate payload does not match its identity.')


def recovery_records(directory):
    """Read the single completion slot and legacy files without modifying either."""
    directory = Path(directory)
    recovery = directory/'recovery.json'
    paths = ([recovery] if recovery.exists() else []) + sorted((directory/'candidates').glob('*.json'))
    for path in paths:
        yield json.loads(path.read_text())
