# Version: 0.29

"""The Versioneer - like a rocketeer, but for versions.
The Versioneer
==============
* like a rocketeer, but for versions!
* https://github.com/python-versioneer/python-versioneer
* Brian Warner
* License: Public Domain (Unlicense)
* Compatible with: Python 3.7, 3.8, 3.9, 3.10, 3.11 and pypy3
* [![Latest Version][pypi-image]][pypi-url]
* [![Build Status][travis-image]][travis-url]
This is a tool for managing a recorded version number in setuptools-based
python projects. The goal is to remove the tedious and error-prone "update
the embedded version string" step from your release process. Making a new
release should be as easy as recording a new tag in your version-control
	@@ -23,38 +27,9 @@
## Quick Install
Versioneer provides two installation modes. The "classic" vendored mode installs
a copy of versioneer into your repository. The experimental build-time dependency mode
is intended to allow you to skip this step and simplify the process of upgrading.
### Vendored mode
* `pip install versioneer` to somewhere in your $PATH
   * A [conda-forge recipe](https://github.com/conda-forge/versioneer-feedstock) is
     available, so you can also use `conda install -c conda-forge versioneer`
* add a `[tool.versioneer]` section to your `pyproject.toml` or a
  `[versioneer]` section to your `setup.cfg` (see [Install](INSTALL.md))
   * Note that you will need to add `tomli; python_version < "3.11"` to your
     build-time dependencies if you use `pyproject.toml`
* run `versioneer install --vendor` in your source tree, commit the results
* verify version information with `python setup.py version`
### Build-time dependency mode
* `pip install versioneer` to somewhere in your $PATH
   * A [conda-forge recipe](https://github.com/conda-forge/versioneer-feedstock) is
     available, so you can also use `conda install -c conda-forge versioneer`
* add a `[tool.versioneer]` section to your `pyproject.toml` or a
  `[versioneer]` section to your `setup.cfg` (see [Install](INSTALL.md))
* add `versioneer` (with `[toml]` extra, if configuring in `pyproject.toml`)
  to the `requires` key of the `build-system` table in `pyproject.toml`:
  ```toml
  [build-system]
  requires = ["setuptools", "versioneer[toml]"]
  build-backend = "setuptools.build_meta"
  ```
* run `versioneer install --no-vendor` in your source tree, commit the results
* verify version information with `python setup.py version`
## Version Identifiers
	@@ -86,7 +61,7 @@
for example `git describe --tags --dirty --always` reports things like
"0.7-1-g574ab98-dirty" to indicate that the checkout is one revision past the
0.7 tag, has a unique revision id of "574ab98", and is "dirty" (it has
uncommitted changes).
The version identifier is used for multiple purposes:
	@@ -191,7 +166,7 @@
Some situations are known to cause problems for Versioneer. This details the
most significant ones. More can be found on Github
[issues page](https://github.com/python-versioneer/python-versioneer/issues).
### Subprojects
	@@ -205,7 +180,7 @@
  `setup.cfg`, and `tox.ini`. Projects like these produce multiple PyPI
  distributions (and upload multiple independently-installable tarballs).
* Source trees whose main purpose is to contain a C library, but which also
  provide bindings to Python (and perhaps other languages) in subdirectories.
Versioneer will look for `.git` in parent directories, and most operations
should get the right version string. However `pip` and `setuptools` have bugs
	@@ -219,9 +194,9 @@
Pip-8.1.1 is known to have this problem, but hopefully it will get fixed in
some later version.
[Bug #38](https://github.com/python-versioneer/python-versioneer/issues/38) is tracking
this issue. The discussion in
[PR #61](https://github.com/python-versioneer/python-versioneer/pull/61) describes the
issue from the Versioneer side in more detail.
[pip PR#3176](https://github.com/pypa/pip/pull/3176) and
[pip PR#3615](https://github.com/pypa/pip/pull/3615) contain work to improve
	@@ -249,20 +224,31 @@
cause egg_info to be rebuilt (including `sdist`, `wheel`, and installing into
a different virtualenv), so this can be surprising.
[Bug #83](https://github.com/python-versioneer/python-versioneer/issues/83) describes
this one, but upgrading to a newer version of setuptools should probably
resolve it.
## Updating Versioneer
To upgrade your project to a new release of Versioneer, do the following:
* install the new Versioneer (`pip install -U versioneer` or equivalent)
* edit `setup.cfg` and `pyproject.toml`, if necessary,
  to include any new configuration settings indicated by the release notes.
  See [UPGRADING](./UPGRADING.md) for details.
* re-run `versioneer install --[no-]vendor` in your source tree, to replace
  `SRC/_version.py`
* commit any changed files
	@@ -279,94 +265,49 @@
direction and include code from all supported VCS systems, reducing the
number of intermediate scripts.
## Similar projects
* [setuptools_scm](https://github.com/pypa/setuptools_scm/) - a non-vendored build-time
  dependency
* [minver](https://github.com/jbweston/miniver) - a lightweight reimplementation of
  versioneer
* [versioningit](https://github.com/jwodder/versioningit) - a PEP 518-based setuptools
  plugin
## License
To make Versioneer easier to embed, all its code is dedicated to the public
domain. The `_version.py` that it creates is also in the public domain.
Specifically, both are released under the "Unlicense", as described in
https://unlicense.org/.
[pypi-image]: https://img.shields.io/pypi/v/versioneer.svg
[pypi-url]: https://pypi.python.org/pypi/versioneer/
[travis-image]:
https://img.shields.io/travis/com/python-versioneer/python-versioneer.svg
[travis-url]: https://travis-ci.com/github/python-versioneer/python-versioneer
"""
# pylint:disable=invalid-name,import-outside-toplevel,missing-function-docstring
# pylint:disable=missing-class-docstring,too-many-branches,too-many-statements
# pylint:disable=raise-missing-from,too-many-lines,too-many-locals,import-error
# pylint:disable=too-few-public-methods,redefined-outer-name,consider-using-with
# pylint:disable=attribute-defined-outside-init,too-many-arguments

import configparser
import errno
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, cast, Dict, List, Optional, Tuple, Union
from typing import NoReturn
import functools

have_tomllib = True
if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomli as tomllib
    except ImportError:
        have_tomllib = False


class VersioneerConfig:
    """Container for Versioneer configuration parameters."""

    VCS: str
    style: str
    tag_prefix: str
    versionfile_source: str
    versionfile_build: Optional[str]
    parentdir_prefix: Optional[str]
    verbose: Optional[bool]


def get_root() -> str:
    """Get the project root directory.
    We require that all commands are run from the project root, i.e. the
    directory that contains setup.py, setup.cfg, and versioneer.py .
    """
    root = os.path.realpath(os.path.abspath(os.getcwd()))
    setup_py = os.path.join(root, "setup.py")
    pyproject_toml = os.path.join(root, "pyproject.toml")
    versioneer_py = os.path.join(root, "versioneer.py")
    if not (
        os.path.exists(setup_py)
        or os.path.exists(pyproject_toml)
        or os.path.exists(versioneer_py)
    ):
        # allow 'python path/to/setup.py COMMAND'
        root = os.path.dirname(os.path.realpath(os.path.abspath(sys.argv[0])))
        setup_py = os.path.join(root, "setup.py")
        pyproject_toml = os.path.join(root, "pyproject.toml")
        versioneer_py = os.path.join(root, "versioneer.py")
    if not (
        os.path.exists(setup_py)
        or os.path.exists(pyproject_toml)
        or os.path.exists(versioneer_py)
    ):
        err = ("Versioneer was unable to run the project root directory. "
               "Versioneer requires setup.py to be executed from "
               "its immediate directory (like 'python setup.py COMMAND'), "
	@@ -380,62 +321,43 @@ def get_root() -> str:
        # module-import table will cache the first one. So we can't use
        # os.path.dirname(__file__), as that will find whichever
        # versioneer.py was first imported, even in later projects.
        my_path = os.path.realpath(os.path.abspath(__file__))
        me_dir = os.path.normcase(os.path.splitext(my_path)[0])
        vsr_dir = os.path.normcase(os.path.splitext(versioneer_py)[0])
        if me_dir != vsr_dir and "VERSIONEER_PEP518" not in globals():
            print("Warning: build in %s is using versioneer.py from %s"
                  % (os.path.dirname(my_path), versioneer_py))
    except NameError:
        pass
    return root


def get_config_from_root(root: str) -> VersioneerConfig:
    """Read the project setup.cfg file to determine Versioneer config."""
    # This might raise OSError (if setup.cfg is missing), or
    # configparser.NoSectionError (if it lacks a [versioneer] section), or
    # configparser.NoOptionError (if it lacks "VCS="). See the docstring at
    # the top of versioneer.py for instructions on writing your setup.cfg .
    root_pth = Path(root)
    pyproject_toml = root_pth / "pyproject.toml"
    setup_cfg = root_pth / "setup.cfg"
    section: Union[Dict[str, Any], configparser.SectionProxy, None] = None
    if pyproject_toml.exists() and have_tomllib:
        try:
            with open(pyproject_toml, 'rb') as fobj:
                pp = tomllib.load(fobj)
            section = pp['tool']['versioneer']
        except (tomllib.TOMLDecodeError, KeyError) as e:
            print(f"Failed to load config from {pyproject_toml}: {e}")
            print("Try to load it from setup.cfg")
    if not section:
        parser = configparser.ConfigParser()
        with open(setup_cfg) as cfg_file:
            parser.read_file(cfg_file)
        parser.get("versioneer", "VCS")  # raise error if missing

        section = parser["versioneer"]

    # `cast`` really shouldn't be used, but its simplest for the
    # common VersioneerConfig users at the moment. We verify against
    # `None` values elsewhere where it matters

    cfg = VersioneerConfig()
    cfg.VCS = section['VCS']
    cfg.style = section.get("style", "")
    cfg.versionfile_source = cast(str, section.get("versionfile_source"))
    cfg.versionfile_build = section.get("versionfile_build")
    cfg.tag_prefix = cast(str, section.get("tag_prefix"))
    if cfg.tag_prefix in ("''", '""', None):
        cfg.tag_prefix = ""
    cfg.parentdir_prefix = section.get("parentdir_prefix")
    if isinstance(section, configparser.SectionProxy):
        # Make sure configparser translates to bool
        cfg.verbose = section.getboolean("verbose")
    else:
        cfg.verbose = section.get("verbose")

    return cfg


	@@ -444,48 +366,37 @@ class NotThisMethod(Exception):


# these dictionaries contain VCS-specific tools
LONG_VERSION_PY: Dict[str, str] = {}
HANDLERS: Dict[str, Dict[str, Callable]] = {}


def register_vcs_handler(vcs: str, method: str) -> Callable:  # decorator
    """Create decorator to mark a method as the handler of a VCS."""
    def decorate(f: Callable) -> Callable:
        """Store f in HANDLERS[vcs][method]."""
        HANDLERS.setdefault(vcs, {})[method] = f
        return f
    return decorate


def run_command(
    commands: List[str],
    args: List[str],
    cwd: Optional[str] = None,
    verbose: bool = False,
    hide_stderr: bool = False,
    env: Optional[Dict[str, str]] = None,
) -> Tuple[Optional[str], Optional[int]]:
    """Call the given command(s)."""
    assert isinstance(commands, list)
    process = None

    popen_kwargs: Dict[str, Any] = {}
    if sys.platform == "win32":
        # This hides the console window if pythonw.exe is used
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        popen_kwargs["startupinfo"] = startupinfo

    for command in commands:
        try:
            dispcmd = str([command] + args)
            # remember shell=False, so use git.cmd on windows, not just git
            process = subprocess.Popen([command] + args, cwd=cwd, env=env,
                                       stdout=subprocess.PIPE,
                                       stderr=(subprocess.PIPE if hide_stderr
                                               else None), **popen_kwargs)
            break
        except OSError as e:
            if e.errno == errno.ENOENT:
                continue
            if verbose:
	@@ -496,25 +407,26 @@ def run_command(
        if verbose:
            print("unable to find command, tried %s" % (commands,))
        return None, None
    stdout = process.communicate()[0].strip().decode()
    if process.returncode != 0:
        if verbose:
            print("unable to run %s (error)" % dispcmd)
            print("stdout was %s" % stdout)
        return None, process.returncode
    return stdout, process.returncode


LONG_VERSION_PY['git'] = r'''
# This file helps to compute a version number in source trees obtained from
# git-archive tarball (such as those provided by githubs download-from-tag
# feature). Distribution tarballs (built by setup.py sdist) and build
# directories (produced by setup.py build) will contain a much shorter file
# that just contains the computed version number.
# This file is released into the public domain.
# Generated by versioneer-0.29
# https://github.com/python-versioneer/python-versioneer
"""Git implementation of _version.py."""
	@@ -523,11 +435,9 @@ def run_command(
import re
import subprocess
import sys
from typing import Any, Callable, Dict, List, Optional, Tuple
import functools
def get_keywords() -> Dict[str, str]:
    """Get the keywords needed to look up the version information."""
    # these strings will be replaced by git during git-archive.
    # setup.py/versioneer.py will grep for the variable names, so they must
	@@ -543,15 +453,8 @@ def get_keywords() -> Dict[str, str]:
class VersioneerConfig:
    """Container for Versioneer configuration parameters."""
    VCS: str
    style: str
    tag_prefix: str
    parentdir_prefix: str
    versionfile_source: str
    verbose: bool
def get_config() -> VersioneerConfig:
    """Create, populate and return the VersioneerConfig() object."""
    # these strings are filled in when 'setup.py versioneer' creates
    # _version.py
	@@ -569,13 +472,13 @@ class NotThisMethod(Exception):
    """Exception raised if a method is not valid for the current scenario."""
LONG_VERSION_PY: Dict[str, str] = {}
HANDLERS: Dict[str, Dict[str, Callable]] = {}
def register_vcs_handler(vcs: str, method: str) -> Callable:  # decorator
    """Create decorator to mark a method as the handler of a VCS."""
    def decorate(f: Callable) -> Callable:
        """Store f in HANDLERS[vcs][method]."""
        if vcs not in HANDLERS:
            HANDLERS[vcs] = {}
	@@ -584,35 +487,22 @@ def decorate(f: Callable) -> Callable:
    return decorate
def run_command(
    commands: List[str],
    args: List[str],
    cwd: Optional[str] = None,
    verbose: bool = False,
    hide_stderr: bool = False,
    env: Optional[Dict[str, str]] = None,
) -> Tuple[Optional[str], Optional[int]]:
    """Call the given command(s)."""
    assert isinstance(commands, list)
    process = None
    popen_kwargs: Dict[str, Any] = {}
    if sys.platform == "win32":
        # This hides the console window if pythonw.exe is used
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        popen_kwargs["startupinfo"] = startupinfo
    for command in commands:
        try:
            dispcmd = str([command] + args)
            # remember shell=False, so use git.cmd on windows, not just git
            process = subprocess.Popen([command] + args, cwd=cwd, env=env,
                                       stdout=subprocess.PIPE,
                                       stderr=(subprocess.PIPE if hide_stderr
                                               else None), **popen_kwargs)
            break
        except OSError as e:
            if e.errno == errno.ENOENT:
                continue
            if verbose:
	@@ -623,20 +513,18 @@ def run_command(
        if verbose:
            print("unable to find command, tried %%s" %% (commands,))
        return None, None
    stdout = process.communicate()[0].strip().decode()
    if process.returncode != 0:
        if verbose:
            print("unable to run %%s (error)" %% dispcmd)
            print("stdout was %%s" %% stdout)
        return None, process.returncode
    return stdout, process.returncode
def versions_from_parentdir(
    parentdir_prefix: str,
    root: str,
    verbose: bool,
) -> Dict[str, Any]:
    """Try to determine the version from the parent directory name.
    Source tarballs conventionally unpack into a directory that includes both
	@@ -645,14 +533,15 @@ def versions_from_parentdir(
    """
    rootdirs = []
    for _ in range(3):
        dirname = os.path.basename(root)
        if dirname.startswith(parentdir_prefix):
            return {"version": dirname[len(parentdir_prefix):],
                    "full-revisionid": None,
                    "dirty": False, "error": None, "date": None}
        rootdirs.append(root)
        root = os.path.dirname(root)  # up a level
    if verbose:
        print("Tried directories %%s but none started with prefix %%s" %%
	@@ -661,48 +550,41 @@ def versions_from_parentdir(
@register_vcs_handler("git", "get_keywords")
def git_get_keywords(versionfile_abs: str) -> Dict[str, str]:
    """Extract version information from the given file."""
    # the code embedded in _version.py can just fetch the value of these
    # keywords. When used from setup.py, we don't want to import _version.py,
    # so we do it with a regexp instead. This function is not used from
    # _version.py.
    keywords: Dict[str, str] = {}
    try:
        with open(versionfile_abs, "r") as fobj:
            for line in fobj:
                if line.strip().startswith("git_refnames ="):
                    mo = re.search(r'=\s*"(.*)"', line)
                    if mo:
                        keywords["refnames"] = mo.group(1)
                if line.strip().startswith("git_full ="):
                    mo = re.search(r'=\s*"(.*)"', line)
                    if mo:
                        keywords["full"] = mo.group(1)
                if line.strip().startswith("git_date ="):
                    mo = re.search(r'=\s*"(.*)"', line)
                    if mo:
                        keywords["date"] = mo.group(1)
    except OSError:
        pass
    return keywords
@register_vcs_handler("git", "keywords")
def git_versions_from_keywords(
    keywords: Dict[str, str],
    tag_prefix: str,
    verbose: bool,
) -> Dict[str, Any]:
    """Get version information from git keywords."""
    if "refnames" not in keywords:
        raise NotThisMethod("Short version file found")
    date = keywords.get("date")
    if date is not None:
        # Use only the last line.  Previous lines may contain GPG signature
        # information.
        date = date.splitlines()[-1]
        # git-2.2.0 added "%%cI", which expands to an ISO-8601 -compliant
        # datestamp. However we prefer "%%ci" (which expands to an "ISO-8601
        # -like" string, which we must then edit to make compliant), because
	@@ -715,11 +597,11 @@ def git_versions_from_keywords(
        if verbose:
            print("keywords are unexpanded, not using")
        raise NotThisMethod("unexpanded keywords, not a git-archive tarball")
    refs = {r.strip() for r in refnames.strip("()").split(",")}
    # starting in git-1.8.3, tags are listed as "tag: foo-1.0" instead of
    # just "foo-1.0". If we see a "tag: " prefix, prefer those.
    TAG = "tag: "
    tags = {r[len(TAG):] for r in refs if r.startswith(TAG)}
    if not tags:
        # Either we're using git < 1.8.3, or there really are no tags. We use
        # a heuristic: assume all version tags have a digit. The old git %%d
	@@ -728,7 +610,7 @@ def git_versions_from_keywords(
        # between branches and tags. By ignoring refnames without digits, we
        # filter out many common branch names like "release" and
        # "stabilization", as well as "HEAD" and "master".
        tags = {r for r in refs if re.search(r'\d', r)}
        if verbose:
            print("discarding '%%s', no digits" %% ",".join(refs - tags))
    if verbose:
	@@ -737,11 +619,6 @@ def git_versions_from_keywords(
        # sorting will prefer e.g. "2.0" over "2.0rc1"
        if ref.startswith(tag_prefix):
            r = ref[len(tag_prefix):]
            # Filter out refs that exactly match prefix or that don't start
            # with a number once the prefix is stripped (mostly a concern
            # when prefix is '')
            if not re.match(r'\d', r):
                continue
            if verbose:
                print("picking %%s" %% r)
            return {"version": r,
	@@ -757,12 +634,7 @@ def git_versions_from_keywords(
@register_vcs_handler("git", "pieces_from_vcs")
def git_pieces_from_vcs(
    tag_prefix: str,
    root: str,
    verbose: bool,
    runner: Callable = run_command
) -> Dict[str, Any]:
    """Get version from 'git describe' in the root of the source tree.
    This only gets called if the git-archive 'subst' keywords were *not*
	@@ -773,73 +645,33 @@ def git_pieces_from_vcs(
    if sys.platform == "win32":
        GITS = ["git.cmd", "git.exe"]
    # GIT_DIR can interfere with correct operation of Versioneer.
    # It may be intended to be passed to the Versioneer-versioned project,
    # but that should not change where we get our version from.
    env = os.environ.copy()
    env.pop("GIT_DIR", None)
    runner = functools.partial(runner, env=env)
    _, rc = runner(GITS, ["rev-parse", "--git-dir"], cwd=root,
                   hide_stderr=not verbose)
    if rc != 0:
        if verbose:
            print("Directory %%s not under git control" %% root)
        raise NotThisMethod("'git rev-parse --git-dir' returned error")
    # if there is a tag matching tag_prefix, this yields TAG-NUM-gHEX[-dirty]
    # if there isn't one, this yields HEX[-dirty] (no NUM)
    describe_out, rc = runner(GITS, [
        "describe", "--tags", "--dirty", "--always", "--long",
        "--match", f"{tag_prefix}[[:digit:]]*"
    ], cwd=root)
    # --long was added in git-1.5.5
    if describe_out is None:
        raise NotThisMethod("'git describe' failed")
    describe_out = describe_out.strip()
    full_out, rc = runner(GITS, ["rev-parse", "HEAD"], cwd=root)
    if full_out is None:
        raise NotThisMethod("'git rev-parse' failed")
    full_out = full_out.strip()
    pieces: Dict[str, Any] = {}
    pieces["long"] = full_out
    pieces["short"] = full_out[:7]  # maybe improved later
    pieces["error"] = None
    branch_name, rc = runner(GITS, ["rev-parse", "--abbrev-ref", "HEAD"],
                             cwd=root)
    # --abbrev-ref was added in git-1.6.3
    if rc != 0 or branch_name is None:
        raise NotThisMethod("'git rev-parse --abbrev-ref' returned error")
    branch_name = branch_name.strip()
    if branch_name == "HEAD":
        # If we aren't exactly on a branch, pick a branch which represents
        # the current commit. If all else fails, we are on a branchless
        # commit.
        branches, rc = runner(GITS, ["branch", "--contains"], cwd=root)
        # --contains was added in git-1.5.4
        if rc != 0 or branches is None:
            raise NotThisMethod("'git branch --contains' returned error")
        branches = branches.split("\n")
        # Remove the first line if we're running detached
        if "(" in branches[0]:
            branches.pop(0)
        # Strip off the leading "* " from the list of branches.
        branches = [branch[2:] for branch in branches]
        if "master" in branches:
            branch_name = "master"
        elif not branches:
            branch_name = None
        else:
            # Pick the first branch that is returned. Good or bad.
            branch_name = branches[0]
    pieces["branch"] = branch_name
    # parse describe_out. It will be like TAG-NUM-gHEX[-dirty] or HEX[-dirty]
    # TAG might have hyphens.
    git_describe = describe_out
	@@ -856,7 +688,7 @@ def git_pieces_from_vcs(
        # TAG-NUM-gHEX
        mo = re.search(r'^(.+)-(\d+)-g([0-9a-f]+)$', git_describe)
        if not mo:
            # unparsable. Maybe git-describe is misbehaving?
            pieces["error"] = ("unable to parse git-describe output: '%%s'"
                               %% describe_out)
            return pieces
	@@ -881,27 +713,26 @@ def git_pieces_from_vcs(
    else:
        # HEX: no tags
        pieces["closest-tag"] = None
        out, rc = runner(GITS, ["rev-list", "HEAD", "--left-right"], cwd=root)
        pieces["distance"] = len(out.split())  # total number of commits
    # commit date: see ISO-8601 comment in git_versions_from_keywords()
    date = runner(GITS, ["show", "-s", "--format=%%ci", "HEAD"], cwd=root)[0].strip()
    # Use only the last line.  Previous lines may contain GPG signature
    # information.
    date = date.splitlines()[-1]
    pieces["date"] = date.strip().replace(" ", "T", 1).replace(" ", "", 1)
    return pieces
def plus_or_dot(pieces: Dict[str, Any]) -> str:
    """Return a + if we don't already have one, else return a ."""
    if "+" in pieces.get("closest-tag", ""):
        return "."
    return "+"
def render_pep440(pieces: Dict[str, Any]) -> str:
    """Build up version string, with post-release "local version identifier".
    Our goal: TAG[+DISTANCE.gHEX[.dirty]] . Note that if you
	@@ -926,71 +757,23 @@ def render_pep440(pieces: Dict[str, Any]) -> str:
    return rendered
def render_pep440_branch(pieces: Dict[str, Any]) -> str:
    """TAG[[.dev0]+DISTANCE.gHEX[.dirty]] .
    The ".dev0" means not master branch. Note that .dev0 sorts backwards
    (a feature branch will appear "older" than the master branch).
    Exceptions:
    1: no tags. 0[.dev0]+untagged.DISTANCE.gHEX[.dirty]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            if pieces["branch"] != "master":
                rendered += ".dev0"
            rendered += plus_or_dot(pieces)
            rendered += "%%d.g%%s" %% (pieces["distance"], pieces["short"])
            if pieces["dirty"]:
                rendered += ".dirty"
    else:
        # exception #1
        rendered = "0"
        if pieces["branch"] != "master":
            rendered += ".dev0"
        rendered += "+untagged.%%d.g%%s" %% (pieces["distance"],
                                          pieces["short"])
        if pieces["dirty"]:
            rendered += ".dirty"
    return rendered
def pep440_split_post(ver: str) -> Tuple[str, Optional[int]]:
    """Split pep440 version string at the post-release segment.
    Returns the release segments before the post-release and the
    post-release version number (or -1 if no post-release segment is present).
    """
    vc = str.split(ver, ".post")
    return vc[0], int(vc[1] or 0) if len(vc) == 2 else None
def render_pep440_pre(pieces: Dict[str, Any]) -> str:
    """TAG[.postN.devDISTANCE] -- No -dirty.
    Exceptions:
    1: no tags. 0.post0.devDISTANCE
    """
    if pieces["closest-tag"]:
        if pieces["distance"]:
            # update the post release segment
            tag_version, post_version = pep440_split_post(pieces["closest-tag"])
            rendered = tag_version
            if post_version is not None:
                rendered += ".post%%d.dev%%d" %% (post_version + 1, pieces["distance"])
            else:
                rendered += ".post0.dev%%d" %% (pieces["distance"])
        else:
            # no commits, use the tag as the version
            rendered = pieces["closest-tag"]
    else:
        # exception #1
        rendered = "0.post0.dev%%d" %% pieces["distance"]
    return rendered
def render_pep440_post(pieces: Dict[str, Any]) -> str:
    """TAG[.postDISTANCE[.dev0]+gHEX] .
    The ".dev0" means dirty. Note that .dev0 sorts backwards
	@@ -1017,41 +800,12 @@ def render_pep440_post(pieces: Dict[str, Any]) -> str:
    return rendered
def render_pep440_post_branch(pieces: Dict[str, Any]) -> str:
    """TAG[.postDISTANCE[.dev0]+gHEX[.dirty]] .
    The ".dev0" means not master branch.
    Exceptions:
    1: no tags. 0.postDISTANCE[.dev0]+gHEX[.dirty]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            rendered += ".post%%d" %% pieces["distance"]
            if pieces["branch"] != "master":
                rendered += ".dev0"
            rendered += plus_or_dot(pieces)
            rendered += "g%%s" %% pieces["short"]
            if pieces["dirty"]:
                rendered += ".dirty"
    else:
        # exception #1
        rendered = "0.post%%d" %% pieces["distance"]
        if pieces["branch"] != "master":
            rendered += ".dev0"
        rendered += "+g%%s" %% pieces["short"]
        if pieces["dirty"]:
            rendered += ".dirty"
    return rendered
def render_pep440_old(pieces: Dict[str, Any]) -> str:
    """TAG[.postDISTANCE[.dev0]] .
    The ".dev0" means dirty.
    Exceptions:
    1: no tags. 0.postDISTANCE[.dev0]
    """
    if pieces["closest-tag"]:
	@@ -1068,7 +822,7 @@ def render_pep440_old(pieces: Dict[str, Any]) -> str:
    return rendered
def render_git_describe(pieces: Dict[str, Any]) -> str:
    """TAG[-DISTANCE-gHEX][-dirty].
    Like 'git describe --tags --dirty --always'.
	@@ -1088,7 +842,7 @@ def render_git_describe(pieces: Dict[str, Any]) -> str:
    return rendered
def render_git_describe_long(pieces: Dict[str, Any]) -> str:
    """TAG-DISTANCE-gHEX[-dirty].
    Like 'git describe --tags --dirty --always -long'.
	@@ -1108,7 +862,7 @@ def render_git_describe_long(pieces: Dict[str, Any]) -> str:
    return rendered
def render(pieces: Dict[str, Any], style: str) -> Dict[str, Any]:
    """Render the given version pieces into the requested style."""
    if pieces["error"]:
        return {"version": "unknown",
	@@ -1122,14 +876,10 @@ def render(pieces: Dict[str, Any], style: str) -> Dict[str, Any]:
    if style == "pep440":
        rendered = render_pep440(pieces)
    elif style == "pep440-branch":
        rendered = render_pep440_branch(pieces)
    elif style == "pep440-pre":
        rendered = render_pep440_pre(pieces)
    elif style == "pep440-post":
        rendered = render_pep440_post(pieces)
    elif style == "pep440-post-branch":
        rendered = render_pep440_post_branch(pieces)
    elif style == "pep440-old":
        rendered = render_pep440_old(pieces)
    elif style == "git-describe":
	@@ -1144,7 +894,7 @@ def render(pieces: Dict[str, Any], style: str) -> Dict[str, Any]:
            "date": pieces.get("date")}
def get_versions() -> Dict[str, Any]:
    """Get version information or return default if unable to do so."""
    # I am in _version.py, which lives at ROOT/VERSIONFILE_SOURCE. If we have
    # __file__, we can work backwards from there to the root. Some
	@@ -1165,7 +915,7 @@ def get_versions() -> Dict[str, Any]:
        # versionfile_source is the relative path from the top of the source
        # tree (where the .git directory might live) to this file. Invert
        # this to find the root from __file__.
        for _ in cfg.versionfile_source.split('/'):
            root = os.path.dirname(root)
    except NameError:
        return {"version": "0+unknown", "full-revisionid": None,
	@@ -1192,48 +942,41 @@ def get_versions() -> Dict[str, Any]:


@register_vcs_handler("git", "get_keywords")
def git_get_keywords(versionfile_abs: str) -> Dict[str, str]:
    """Extract version information from the given file."""
    # the code embedded in _version.py can just fetch the value of these
    # keywords. When used from setup.py, we don't want to import _version.py,
    # so we do it with a regexp instead. This function is not used from
    # _version.py.
    keywords: Dict[str, str] = {}
    try:
        with open(versionfile_abs, "r") as fobj:
            for line in fobj:
                if line.strip().startswith("git_refnames ="):
                    mo = re.search(r'=\s*"(.*)"', line)
                    if mo:
                        keywords["refnames"] = mo.group(1)
                if line.strip().startswith("git_full ="):
                    mo = re.search(r'=\s*"(.*)"', line)
                    if mo:
                        keywords["full"] = mo.group(1)
                if line.strip().startswith("git_date ="):
                    mo = re.search(r'=\s*"(.*)"', line)
                    if mo:
                        keywords["date"] = mo.group(1)
    except OSError:
        pass
    return keywords


@register_vcs_handler("git", "keywords")
def git_versions_from_keywords(
    keywords: Dict[str, str],
    tag_prefix: str,
    verbose: bool,
) -> Dict[str, Any]:
    """Get version information from git keywords."""
    if "refnames" not in keywords:
        raise NotThisMethod("Short version file found")
    date = keywords.get("date")
    if date is not None:
        # Use only the last line.  Previous lines may contain GPG signature
        # information.
        date = date.splitlines()[-1]

        # git-2.2.0 added "%cI", which expands to an ISO-8601 -compliant
        # datestamp. However we prefer "%ci" (which expands to an "ISO-8601
        # -like" string, which we must then edit to make compliant), because
	@@ -1246,11 +989,11 @@ def git_versions_from_keywords(
        if verbose:
            print("keywords are unexpanded, not using")
        raise NotThisMethod("unexpanded keywords, not a git-archive tarball")
    refs = {r.strip() for r in refnames.strip("()").split(",")}
    # starting in git-1.8.3, tags are listed as "tag: foo-1.0" instead of
    # just "foo-1.0". If we see a "tag: " prefix, prefer those.
    TAG = "tag: "
    tags = {r[len(TAG):] for r in refs if r.startswith(TAG)}
    if not tags:
        # Either we're using git < 1.8.3, or there really are no tags. We use
        # a heuristic: assume all version tags have a digit. The old git %d
	@@ -1259,7 +1002,7 @@ def git_versions_from_keywords(
        # between branches and tags. By ignoring refnames without digits, we
        # filter out many common branch names like "release" and
        # "stabilization", as well as "HEAD" and "master".
        tags = {r for r in refs if re.search(r'\d', r)}
        if verbose:
            print("discarding '%s', no digits" % ",".join(refs - tags))
    if verbose:
	@@ -1268,11 +1011,6 @@ def git_versions_from_keywords(
        # sorting will prefer e.g. "2.0" over "2.0rc1"
        if ref.startswith(tag_prefix):
            r = ref[len(tag_prefix):]
            # Filter out refs that exactly match prefix or that don't start
            # with a number once the prefix is stripped (mostly a concern
            # when prefix is '')
            if not re.match(r'\d', r):
                continue
            if verbose:
                print("picking %s" % r)
            return {"version": r,
	@@ -1288,12 +1026,7 @@ def git_versions_from_keywords(


@register_vcs_handler("git", "pieces_from_vcs")
def git_pieces_from_vcs(
    tag_prefix: str,
    root: str,
    verbose: bool,
    runner: Callable = run_command
) -> Dict[str, Any]:
    """Get version from 'git describe' in the root of the source tree.
    This only gets called if the git-archive 'subst' keywords were *not*
	@@ -1304,73 +1037,33 @@ def git_pieces_from_vcs(
    if sys.platform == "win32":
        GITS = ["git.cmd", "git.exe"]

    # GIT_DIR can interfere with correct operation of Versioneer.
    # It may be intended to be passed to the Versioneer-versioned project,
    # but that should not change where we get our version from.
    env = os.environ.copy()
    env.pop("GIT_DIR", None)
    runner = functools.partial(runner, env=env)

    _, rc = runner(GITS, ["rev-parse", "--git-dir"], cwd=root,
                   hide_stderr=not verbose)
    if rc != 0:
        if verbose:
            print("Directory %s not under git control" % root)
        raise NotThisMethod("'git rev-parse --git-dir' returned error")

    # if there is a tag matching tag_prefix, this yields TAG-NUM-gHEX[-dirty]
    # if there isn't one, this yields HEX[-dirty] (no NUM)
    describe_out, rc = runner(GITS, [
        "describe", "--tags", "--dirty", "--always", "--long",
        "--match", f"{tag_prefix}[[:digit:]]*"
    ], cwd=root)
    # --long was added in git-1.5.5
    if describe_out is None:
        raise NotThisMethod("'git describe' failed")
    describe_out = describe_out.strip()
    full_out, rc = runner(GITS, ["rev-parse", "HEAD"], cwd=root)
    if full_out is None:
        raise NotThisMethod("'git rev-parse' failed")
    full_out = full_out.strip()

    pieces: Dict[str, Any] = {}
    pieces["long"] = full_out
    pieces["short"] = full_out[:7]  # maybe improved later
    pieces["error"] = None

    branch_name, rc = runner(GITS, ["rev-parse", "--abbrev-ref", "HEAD"],
                             cwd=root)
    # --abbrev-ref was added in git-1.6.3
    if rc != 0 or branch_name is None:
        raise NotThisMethod("'git rev-parse --abbrev-ref' returned error")
    branch_name = branch_name.strip()

    if branch_name == "HEAD":
        # If we aren't exactly on a branch, pick a branch which represents
        # the current commit. If all else fails, we are on a branchless
        # commit.
        branches, rc = runner(GITS, ["branch", "--contains"], cwd=root)
        # --contains was added in git-1.5.4
        if rc != 0 or branches is None:
            raise NotThisMethod("'git branch --contains' returned error")
        branches = branches.split("\n")

        # Remove the first line if we're running detached
        if "(" in branches[0]:
            branches.pop(0)

        # Strip off the leading "* " from the list of branches.
        branches = [branch[2:] for branch in branches]
        if "master" in branches:
            branch_name = "master"
        elif not branches:
            branch_name = None
        else:
            # Pick the first branch that is returned. Good or bad.
            branch_name = branches[0]

    pieces["branch"] = branch_name

    # parse describe_out. It will be like TAG-NUM-gHEX[-dirty] or HEX[-dirty]
    # TAG might have hyphens.
    git_describe = describe_out
	@@ -1387,7 +1080,7 @@ def git_pieces_from_vcs(
        # TAG-NUM-gHEX
        mo = re.search(r'^(.+)-(\d+)-g([0-9a-f]+)$', git_describe)
        if not mo:
            # unparsable. Maybe git-describe is misbehaving?
            pieces["error"] = ("unable to parse git-describe output: '%s'"
                               % describe_out)
            return pieces
	@@ -1412,20 +1105,19 @@ def git_pieces_from_vcs(
    else:
        # HEX: no tags
        pieces["closest-tag"] = None
        out, rc = runner(GITS, ["rev-list", "HEAD", "--left-right"], cwd=root)
        pieces["distance"] = len(out.split())  # total number of commits

    # commit date: see ISO-8601 comment in git_versions_from_keywords()
    date = runner(GITS, ["show", "-s", "--format=%ci", "HEAD"], cwd=root)[0].strip()
    # Use only the last line.  Previous lines may contain GPG signature
    # information.
    date = date.splitlines()[-1]
    pieces["date"] = date.strip().replace(" ", "T", 1).replace(" ", "", 1)

    return pieces


def do_vcs_install(versionfile_source: str, ipy: Optional[str]) -> None:
    """Git-specific installation logic for Versioneer.
    For Git, this means creating/changing .gitattributes to mark _version.py
	@@ -1434,40 +1126,36 @@ def do_vcs_install(versionfile_source: str, ipy: Optional[str]) -> None:
    GITS = ["git"]
    if sys.platform == "win32":
        GITS = ["git.cmd", "git.exe"]
    files = [versionfile_source]
    if ipy:
        files.append(ipy)
    if "VERSIONEER_PEP518" not in globals():
        try:
            my_path = __file__
            if my_path.endswith((".pyc", ".pyo")):
                my_path = os.path.splitext(my_path)[0] + ".py"
            versioneer_file = os.path.relpath(my_path)
        except NameError:
            versioneer_file = "versioneer.py"
        files.append(versioneer_file)
    present = False
    try:
        with open(".gitattributes", "r") as fobj:
            for line in fobj:
                if line.strip().startswith(versionfile_source):
                    if "export-subst" in line.strip().split()[1:]:
                        present = True
                        break
    except OSError:
        pass
    if not present:
        with open(".gitattributes", "a+") as fobj:
            fobj.write(f"{versionfile_source} export-subst\n")
        files.append(".gitattributes")
    run_command(GITS, ["add", "--"] + files)


def versions_from_parentdir(
    parentdir_prefix: str,
    root: str,
    verbose: bool,
) -> Dict[str, Any]:
    """Try to determine the version from the parent directory name.
    Source tarballs conventionally unpack into a directory that includes both
	@@ -1476,14 +1164,15 @@ def versions_from_parentdir(
    """
    rootdirs = []

    for _ in range(3):
        dirname = os.path.basename(root)
        if dirname.startswith(parentdir_prefix):
            return {"version": dirname[len(parentdir_prefix):],
                    "full-revisionid": None,
                    "dirty": False, "error": None, "date": None}
        rootdirs.append(root)
        root = os.path.dirname(root)  # up a level

    if verbose:
        print("Tried directories %s but none started with prefix %s" %
	@@ -1492,7 +1181,7 @@ def versions_from_parentdir(


SHORT_VERSION_PY = """
# This file was generated by 'versioneer.py' (0.29) from
# revision-control system data, or from the parent directory name of an
# unpacked source archive. Distribution tarballs contain a pre-generated copy
# of this file.
	@@ -1509,12 +1198,12 @@ def get_versions():
"""


def versions_from_file(filename: str) -> Dict[str, Any]:
    """Try to determine the version from _version.py if present."""
    try:
        with open(filename) as f:
            contents = f.read()
    except OSError:
        raise NotThisMethod("unable to read _version.py")
    mo = re.search(r"version_json = '''\n(.*)'''  # END VERSION_JSON",
                   contents, re.M | re.S)
	@@ -1526,8 +1215,9 @@ def versions_from_file(filename: str) -> Dict[str, Any]:
    return json.loads(mo.group(1))


def write_to_version_file(filename: str, versions: Dict[str, Any]) -> None:
    """Write the given version number to the given _version.py file."""
    contents = json.dumps(versions, sort_keys=True,
                          indent=1, separators=(",", ": "))
    with open(filename, "w") as f:
	@@ -1536,14 +1226,14 @@ def write_to_version_file(filename: str, versions: Dict[str, Any]) -> None:
    print("set %s to '%s'" % (filename, versions["version"]))


def plus_or_dot(pieces: Dict[str, Any]) -> str:
    """Return a + if we don't already have one, else return a ."""
    if "+" in pieces.get("closest-tag", ""):
        return "."
    return "+"


def render_pep440(pieces: Dict[str, Any]) -> str:
    """Build up version string, with post-release "local version identifier".
    Our goal: TAG[+DISTANCE.gHEX[.dirty]] . Note that if you
	@@ -1568,71 +1258,23 @@ def render_pep440(pieces: Dict[str, Any]) -> str:
    return rendered


def render_pep440_branch(pieces: Dict[str, Any]) -> str:
    """TAG[[.dev0]+DISTANCE.gHEX[.dirty]] .
    The ".dev0" means not master branch. Note that .dev0 sorts backwards
    (a feature branch will appear "older" than the master branch).
    Exceptions:
    1: no tags. 0[.dev0]+untagged.DISTANCE.gHEX[.dirty]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            if pieces["branch"] != "master":
                rendered += ".dev0"
            rendered += plus_or_dot(pieces)
            rendered += "%d.g%s" % (pieces["distance"], pieces["short"])
            if pieces["dirty"]:
                rendered += ".dirty"
    else:
        # exception #1
        rendered = "0"
        if pieces["branch"] != "master":
            rendered += ".dev0"
        rendered += "+untagged.%d.g%s" % (pieces["distance"],
                                          pieces["short"])
        if pieces["dirty"]:
            rendered += ".dirty"
    return rendered


def pep440_split_post(ver: str) -> Tuple[str, Optional[int]]:
    """Split pep440 version string at the post-release segment.
    Returns the release segments before the post-release and the
    post-release version number (or -1 if no post-release segment is present).
    """
    vc = str.split(ver, ".post")
    return vc[0], int(vc[1] or 0) if len(vc) == 2 else None


def render_pep440_pre(pieces: Dict[str, Any]) -> str:
    """TAG[.postN.devDISTANCE] -- No -dirty.
    Exceptions:
    1: no tags. 0.post0.devDISTANCE
    """
    if pieces["closest-tag"]:
        if pieces["distance"]:
            # update the post release segment
            tag_version, post_version = pep440_split_post(pieces["closest-tag"])
            rendered = tag_version
            if post_version is not None:
                rendered += ".post%d.dev%d" % (post_version + 1, pieces["distance"])
            else:
                rendered += ".post0.dev%d" % (pieces["distance"])
        else:
            # no commits, use the tag as the version
            rendered = pieces["closest-tag"]
    else:
        # exception #1
        rendered = "0.post0.dev%d" % pieces["distance"]
    return rendered


def render_pep440_post(pieces: Dict[str, Any]) -> str:
    """TAG[.postDISTANCE[.dev0]+gHEX] .
    The ".dev0" means dirty. Note that .dev0 sorts backwards
	@@ -1659,41 +1301,12 @@ def render_pep440_post(pieces: Dict[str, Any]) -> str:
    return rendered


def render_pep440_post_branch(pieces: Dict[str, Any]) -> str:
    """TAG[.postDISTANCE[.dev0]+gHEX[.dirty]] .
    The ".dev0" means not master branch.
    Exceptions:
    1: no tags. 0.postDISTANCE[.dev0]+gHEX[.dirty]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            rendered += ".post%d" % pieces["distance"]
            if pieces["branch"] != "master":
                rendered += ".dev0"
            rendered += plus_or_dot(pieces)
            rendered += "g%s" % pieces["short"]
            if pieces["dirty"]:
                rendered += ".dirty"
    else:
        # exception #1
        rendered = "0.post%d" % pieces["distance"]
        if pieces["branch"] != "master":
            rendered += ".dev0"
        rendered += "+g%s" % pieces["short"]
        if pieces["dirty"]:
            rendered += ".dirty"
    return rendered


def render_pep440_old(pieces: Dict[str, Any]) -> str:
    """TAG[.postDISTANCE[.dev0]] .
    The ".dev0" means dirty.
    Exceptions:
    1: no tags. 0.postDISTANCE[.dev0]
    """
    if pieces["closest-tag"]:
	@@ -1710,7 +1323,7 @@ def render_pep440_old(pieces: Dict[str, Any]) -> str:
    return rendered


def render_git_describe(pieces: Dict[str, Any]) -> str:
    """TAG[-DISTANCE-gHEX][-dirty].
    Like 'git describe --tags --dirty --always'.
	@@ -1730,7 +1343,7 @@ def render_git_describe(pieces: Dict[str, Any]) -> str:
    return rendered


def render_git_describe_long(pieces: Dict[str, Any]) -> str:
    """TAG-DISTANCE-gHEX[-dirty].
    Like 'git describe --tags --dirty --always -long'.
	@@ -1750,7 +1363,7 @@ def render_git_describe_long(pieces: Dict[str, Any]) -> str:
    return rendered


def render(pieces: Dict[str, Any], style: str) -> Dict[str, Any]:
    """Render the given version pieces into the requested style."""
    if pieces["error"]:
        return {"version": "unknown",
	@@ -1764,14 +1377,10 @@ def render(pieces: Dict[str, Any], style: str) -> Dict[str, Any]:

    if style == "pep440":
        rendered = render_pep440(pieces)
    elif style == "pep440-branch":
        rendered = render_pep440_branch(pieces)
    elif style == "pep440-pre":
        rendered = render_pep440_pre(pieces)
    elif style == "pep440-post":
        rendered = render_pep440_post(pieces)
    elif style == "pep440-post-branch":
        rendered = render_pep440_post_branch(pieces)
    elif style == "pep440-old":
        rendered = render_pep440_old(pieces)
    elif style == "git-describe":
	@@ -1790,7 +1399,7 @@ class VersioneerBadRootError(Exception):
    """The project root directory is unknown or missing key files."""


def get_versions(verbose: bool = False) -> Dict[str, Any]:
    """Get the project version from whatever source is available.
    Returns dict with two keys: 'version' and 'full'.
	@@ -1805,7 +1414,7 @@ def get_versions(verbose: bool = False) -> Dict[str, Any]:
    assert cfg.VCS is not None, "please set [versioneer]VCS= in setup.cfg"
    handlers = HANDLERS.get(cfg.VCS)
    assert handlers, "unrecognized VCS '%s'" % cfg.VCS
    verbose = verbose or bool(cfg.verbose)  # `bool()` used to avoid `None`
    assert cfg.versionfile_source is not None, \
        "please set versioneer.versionfile_source"
    assert cfg.tag_prefix is not None, "please set versioneer.tag_prefix"
	@@ -1866,17 +1475,13 @@ def get_versions(verbose: bool = False) -> Dict[str, Any]:
            "date": None}


def get_version() -> str:
    """Get the short version string for this project."""
    return get_versions()["version"]


def get_cmdclass(cmdclass: Optional[Dict[str, Any]] = None):
    """Get the custom setuptools subclasses used by Versioneer.
    If the package uses a different cmdclass (e.g. one from numpy), it
    should be provide as an argument.
    """
    if "versioneer" in sys.modules:
        del sys.modules["versioneer"]
        # this fixes the "python setup.py develop" case (also 'install' and
	@@ -1890,25 +1495,25 @@ def get_cmdclass(cmdclass: Optional[Dict[str, Any]] = None):
        # parent is protected against the child's "import versioneer". By
        # removing ourselves from sys.modules here, before the child build
        # happens, we protect the child from the parent's versioneer too.
        # Also see https://github.com/python-versioneer/python-versioneer/issues/52

    cmds = {} if cmdclass is None else cmdclass.copy()

    # we add "version" to setuptools
    from setuptools import Command

    class cmd_version(Command):
        description = "report generated version string"
        user_options: List[Tuple[str, str, str]] = []
        boolean_options: List[str] = []

        def initialize_options(self) -> None:
            pass

        def finalize_options(self) -> None:
            pass

        def run(self) -> None:
            vers = get_versions(verbose=True)
            print("Version: %s" % vers["version"])
            print(" full-revisionid: %s" % vers.get("full-revisionid"))
	@@ -1918,7 +1523,7 @@ def run(self) -> None:
                print(" error: %s" % vers["error"])
    cmds["version"] = cmd_version

    # we override "build_py" in setuptools
    #
    # most invocation pathways end up running build_py:
    #  distutils/build -> build_py
	@@ -1933,25 +1538,18 @@ def run(self) -> None:
    #   then does setup.py bdist_wheel, or sometimes setup.py install
    #  setup.py egg_info -> ?

    # pip install -e . and setuptool/editable_wheel will invoke build_py
    # but the build_py command is not expected to copy any files.

    # we override different "build_py" commands for both environments
    if 'build_py' in cmds:
        _build_py: Any = cmds['build_py']
    else:
        from setuptools.command.build_py import build_py as _build_py

    class cmd_build_py(_build_py):
        def run(self) -> None:
            root = get_root()
            cfg = get_config_from_root(root)
            versions = get_versions()
            _build_py.run(self)
            if getattr(self, "editable_mode", False):
                # During editable installs `.py` and data files are
                # not copied to build_lib
                return
            # now locate _version.py in the new build/ directory and replace
            # it with an updated value
            if cfg.versionfile_build:
	@@ -1961,40 +1559,8 @@ def run(self) -> None:
                write_to_version_file(target_versionfile, versions)
    cmds["build_py"] = cmd_build_py

    if 'build_ext' in cmds:
        _build_ext: Any = cmds['build_ext']
    else:
        from setuptools.command.build_ext import build_ext as _build_ext

    class cmd_build_ext(_build_ext):
        def run(self) -> None:
            root = get_root()
            cfg = get_config_from_root(root)
            versions = get_versions()
            _build_ext.run(self)
            if self.inplace:
                # build_ext --inplace will only build extensions in
                # build/lib<..> dir with no _version.py to write to.
                # As in place builds will already have a _version.py
                # in the module dir, we do not need to write one.
                return
            # now locate _version.py in the new build/ directory and replace
            # it with an updated value
            if not cfg.versionfile_build:
                return
            target_versionfile = os.path.join(self.build_lib,
                                              cfg.versionfile_build)
            if not os.path.exists(target_versionfile):
                print(f"Warning: {target_versionfile} does not exist, skipping "
                      "version update. This can happen if you are running build_ext "
                      "without first running build_py.")
                return
            print("UPDATING %s" % target_versionfile)
            write_to_version_file(target_versionfile, versions)
    cmds["build_ext"] = cmd_build_ext

    if "cx_Freeze" in sys.modules:  # cx_freeze enabled?
        from cx_Freeze.dist import build_exe as _build_exe  # type: ignore
        # nczeczulin reports that py2exe won't like the pep440-style string
        # as FILEVERSION, but it can be used for PRODUCTVERSION, e.g.
        # setup(console=[{
	@@ -2003,7 +1569,7 @@ def run(self) -> None:
        #   ...

        class cmd_build_exe(_build_exe):
            def run(self) -> None:
                root = get_root()
                cfg = get_config_from_root(root)
                versions = get_versions()
	@@ -2027,12 +1593,12 @@ def run(self) -> None:

    if 'py2exe' in sys.modules:  # py2exe enabled?
        try:
            from py2exe.setuptools_buildexe import py2exe as _py2exe  # type: ignore
        except ImportError:
            from py2exe.distutils_buildexe import py2exe as _py2exe  # type: ignore

        class cmd_py2exe(_py2exe):
            def run(self) -> None:
                root = get_root()
                cfg = get_config_from_root(root)
                versions = get_versions()
	@@ -2053,59 +1619,22 @@ def run(self) -> None:
                             })
        cmds["py2exe"] = cmd_py2exe

    # sdist farms its file list building out to egg_info
    if 'egg_info' in cmds:
        _egg_info: Any = cmds['egg_info']
    else:
        from setuptools.command.egg_info import egg_info as _egg_info

    class cmd_egg_info(_egg_info):
        def find_sources(self) -> None:
            # egg_info.find_sources builds the manifest list and writes it
            # in one shot
            super().find_sources()

            # Modify the filelist and normalize it
            root = get_root()
            cfg = get_config_from_root(root)
            self.filelist.append('versioneer.py')
            if cfg.versionfile_source:
                # There are rare cases where versionfile_source might not be
                # included by default, so we must be explicit
                self.filelist.append(cfg.versionfile_source)
            self.filelist.sort()
            self.filelist.remove_duplicates()

            # The write method is hidden in the manifest_maker instance that
            # generated the filelist and was thrown away
            # We will instead replicate their final normalization (to unicode,
            # and POSIX-style paths)
            from setuptools import unicode_utils
            normalized = [unicode_utils.filesys_decode(f).replace(os.sep, '/')
                          for f in self.filelist.files]

            manifest_filename = os.path.join(self.egg_info, 'SOURCES.txt')
            with open(manifest_filename, 'w') as fobj:
                fobj.write('\n'.join(normalized))

    cmds['egg_info'] = cmd_egg_info

    # we override different "sdist" commands for both environments
    if 'sdist' in cmds:
        _sdist: Any = cmds['sdist']
    else:
        from setuptools.command.sdist import sdist as _sdist

    class cmd_sdist(_sdist):
        def run(self) -> None:
            versions = get_versions()
            self._versioneer_generated_versions = versions
            # unless we update this, the command will keep using the old
            # version
            self.distribution.metadata.version = versions["version"]
            return _sdist.run(self)

        def make_release_tree(self, base_dir: str, files: List[str]) -> None:
            root = get_root()
            cfg = get_config_from_root(root)
            _sdist.make_release_tree(self, base_dir, files)
	@@ -2158,26 +1687,21 @@ def make_release_tree(self, base_dir: str, files: List[str]) -> None:
"""

OLD_SNIPPET = """
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
"""

INIT_PY_SNIPPET = """
from . import {0}
__version__ = {0}.get_versions()['version']
"""


def do_setup() -> int:
    """Do main VCS-independent setup function for installing Versioneer."""
    root = get_root()
    try:
        cfg = get_config_from_root(root)
    except (OSError, configparser.NoSectionError,
            configparser.NoOptionError) as e:
        if isinstance(e, (OSError, configparser.NoSectionError)):
            print("Adding sample versioneer config to setup.cfg",
                  file=sys.stderr)
            with open(os.path.join(root, "setup.cfg"), "a") as f:
	@@ -2197,37 +1721,62 @@ def do_setup() -> int:

    ipy = os.path.join(os.path.dirname(cfg.versionfile_source),
                       "__init__.py")
    maybe_ipy: Optional[str] = ipy
    if os.path.exists(ipy):
        try:
            with open(ipy, "r") as f:
                old = f.read()
        except OSError:
            old = ""
        module = os.path.splitext(os.path.basename(cfg.versionfile_source))[0]
        snippet = INIT_PY_SNIPPET.format(module)
        if OLD_SNIPPET in old:
            print(" replacing boilerplate in %s" % ipy)
            with open(ipy, "w") as f:
                f.write(old.replace(OLD_SNIPPET, snippet))
        elif snippet not in old:
            print(" appending to %s" % ipy)
            with open(ipy, "a") as f:
                f.write(snippet)
        else:
            print(" %s unmodified" % ipy)
    else:
        print(" %s doesn't exist, ok" % ipy)
        maybe_ipy = None

    # Make VCS-specific changes. For git, this means creating/changing
    # .gitattributes to mark _version.py for export-subst keyword
    # substitution.
    do_vcs_install(cfg.versionfile_source, maybe_ipy)
    return 0


def scan_setup_py() -> int:
    """Validate the contents of setup.py against Versioneer's expectations."""
    found = set()
    setters = False
	@@ -2264,14 +1813,10 @@ def scan_setup_py() -> int:
    return errors


def setup_command() -> NoReturn:
    """Set up Versioneer and exit with appropriate error code."""
    errors = do_setup()
    errors += scan_setup_py()
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "setup":
        setup_command()
