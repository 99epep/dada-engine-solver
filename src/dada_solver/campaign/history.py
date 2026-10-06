"""Durable files with an append-only evaluation journal and recoverable snapshots."""
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import fcntl
import gzip
import zlib
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


def journal_path(directory):
    """Use the compressed journal; reject unsupported persistence layouts."""
    directory = Path(directory)
    if (directory / 'history.jsonl').exists():
        raise ValueError('Uncompressed campaign journals are unsupported.')
    if any((directory / 'candidates').glob('*.json')):
        raise ValueError('Per-candidate recovery files are unsupported.')
    return directory / 'history.jsonl.gz'


def read_journal(path):
    """Read a snapshot without repairing it; return records, bytes and torn-tail offset.

    Each gzip member is one durable JSON record. CRC errors are corruption, not
    interrupted writes. Only an incomplete final member may be recovered.
    """
    path=Path(path)
    if path.suffix != '.gz':
        raise ValueError('Campaign journals must use the compressed format.')
    data=path.read_bytes() if path.exists() else b''
    records=[];offset=0
    while offset<len(data):
        start=offset;decoder=zlib.decompressobj(31);chunks=[]
        while offset<len(data) and not decoder.eof:
            chunk=data[offset:offset+65536]
            try: chunks.append(decoder.decompress(chunk))
            except zlib.error as error:
                raise ValueError(f'Corrupt compressed history member at byte {start}.') from error
            offset+=len(chunk)-len(decoder.unused_data)
        if not decoder.eof: return records,data,start
        payload=b''.join(chunks)
        try:
            if not payload.endswith(b'\n'): raise ValueError('Missing record terminator.')
            record=json.loads(payload)
        except (ValueError,UnicodeDecodeError) as error:
            raise ValueError(f'Invalid JSON in complete history member at byte {start}.') from error
        records.append(record)
    return records,data,None


class CampaignHistory:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.path = journal_path(self.directory)

    @contextmanager
    def locked(self):
        with (self.directory/'.writer.lock').open('a') as stream:
            try: fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError('Another writer is already running this campaign.') from error
            try: yield
            finally: fcntl.flock(stream, fcntl.LOCK_UN)

    def load(self):
        records,data,torn_offset=read_journal(self.path)
        if torn_offset is not None:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
            (self.directory/f'history_torn_tail_{stamp}.bin').write_bytes(data[torn_offset:])
            with self.path.open('r+b') as stream:
                stream.truncate(torn_offset); stream.flush(); os.fsync(stream.fileno())
        for record in records: verify_record(record)
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
        if [r['evaluation_number'] for r in records] != list(range(len(records))):
            raise ValueError('History evaluation numbers are not contiguous.')
        return records

    def _append(self, record):
        payload=(canonical_json(record)+'\n').encode('utf-8')
        if self.path.suffix=='.gz':
            payload=gzip.compress(payload,compresslevel=6,mtime=0)
        with self.path.open('ab') as stream:
            stream.write(payload); stream.flush(); os.fsync(stream.fileno())

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
    """Read the single durable completion slot without modifying it."""
    directory = Path(directory)
    recovery = directory/'recovery.json'
    paths = [recovery] if recovery.exists() else []
    for path in paths:
        try: text = path.read_text()
        except FileNotFoundError:
            if path == recovery: continue  # A live writer may have acknowledged it.
            raise
        yield json.loads(text)
