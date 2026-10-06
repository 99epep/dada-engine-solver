"""Distribution boundaries are verified against archives built in isolation."""

from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[1]
APPLICATION_DIRECTORY = 'exam' + 'ples'
EXCLUDED_DIRECTORIES = {APPLICATION_DIRECTORY, 'outputs'}


@pytest.fixture(scope='module')
def distributions(tmp_path_factory):
    stage = tmp_path_factory.mktemp('distribution')
    for name in ('src', 'docs', 'tests'):
        shutil.copytree(ROOT / name, stage / name,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ('LICENSE', 'README.md', 'pyproject.toml', 'MANIFEST.in'):
        shutil.copy2(ROOT / name, stage / name)
    # Exclusions must work even with files present and a stale manifest cache.
    cached_sources = stage / 'src/dada_engine_solver.egg-info/SOURCES.txt'
    cached_sources.parent.mkdir(exist_ok=True)
    with cached_sources.open('a') as stream:
        for directory in EXCLUDED_DIRECTORIES:
            target = stage / directory / 'distribution_sentinel.py'
            target.parent.mkdir(exist_ok=True)
            target.write_text('# Not a distribution resource.\n')
            stream.write(f'\n{target.relative_to(stage).as_posix()}\n')
    result = subprocess.run(
        [sys.executable, '-c',
         'from setuptools.build_meta import build_sdist; build_sdist("dist")'],
        cwd=stage, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    sdist = next((stage / 'dist').glob('*.tar.gz'))
    unpacked = stage / 'unpacked'
    with tarfile.open(sdist) as archive:
        archive.extractall(unpacked, filter='data')
    source = next(unpacked.iterdir())
    result = subprocess.run(
        [sys.executable, '-c',
         'from setuptools.build_meta import build_wheel; build_wheel("dist")'],
        cwd=source, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return sdist, next((source / 'dist').glob('*.whl'))


def test_distribution_excludes_application_and_output_directories(distributions):
    sdist, wheel = distributions
    with tarfile.open(sdist) as archive:
        source_paths = [Path(name).parts[1:] for name in archive.getnames()]
    with zipfile.ZipFile(wheel) as archive:
        wheel_paths = [Path(name).parts for name in archive.namelist()]
    assert not any(set(parts) & EXCLUDED_DIRECTORIES
                   for parts in source_paths + wheel_paths)


def test_distribution_preserves_runtime_resources_and_doty_sources(distributions):
    sdist, wheel = distributions
    package = ROOT / 'src/dada_solver'
    resources = [path.relative_to(ROOT / 'src').as_posix()
                 for path in package.rglob('*')
                 if path.is_file() and (path.suffix == '.py'
                    or (path.parent == package / 'research/data'
                        and path.suffix in {'.json', '.toml'})
                    or (path.parent == package / 'research' and path.suffix == '.js'))]
    with tarfile.open(sdist) as archive:
        members = {name.split('/', 1)[-1]: name for name in archive.getnames()}
        assert all('src/' + name in members for name in resources)
        for species in ('nitrogen', 'helium'):
            name = f'tests/data/doty_1991_{species}_reference.csv'
            assert archive.extractfile(members[name]).read() == (ROOT / name).read_bytes()
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        assert set(resources) <= names
        assert not any(name.startswith('tests/') for name in names)


def test_wheel_can_be_installed_and_load_entry_points(distributions, tmp_path):
    _, wheel = distributions
    target = tmp_path / 'installed'
    install = subprocess.run(
        [sys.executable, '-m', 'pip', 'install', '--no-deps', '--no-index',
         '--target', str(target), str(wheel)], capture_output=True, text=True,
    )
    assert install.returncode == 0, install.stdout + install.stderr
    result = subprocess.run(
        [sys.executable, '-I', '-c',
         'import sys; sys.path.insert(0, sys.argv[1]); '
         'from importlib.metadata import distribution; '
         'from importlib.resources import files; '
         'assert files("dada_solver.research").joinpath("data/machine_basis_v2.json").is_file(); '
         '[entry.load() for entry in distribution("dada-engine-solver").entry_points]',
         str(target)], cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
