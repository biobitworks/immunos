# IMMUNOS Packaging Standards

## Overview

IMMUNOS must comply with packaging standards for major software distribution channels. This document outlines requirements for each platform.

---

## PyPI (pip)

### Required Files

```
immunos/
├── pyproject.toml          # Modern Python packaging (PEP 517/518)
├── setup.py                # Legacy compatibility (optional)
├── setup.cfg               # Configuration (optional, use pyproject.toml)
├── MANIFEST.in             # Include non-Python files
├── LICENSE                 # License file (required)
├── README.md               # Package description
├── CHANGELOG.md            # Version history
├── src/
│   └── immunos/
│       ├── __init__.py     # Package init with __version__
│       ├── cli.py          # Entry points
│       └── ...
└── tests/
    └── ...
```

### pyproject.toml (PEP 621)

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "immunos"
version = "0.1.0"
description = "AI-assisted research integrity system using immune system metaphors"
readme = "README.md"
license = {text = "MIT"}
requires-python = ">=3.9"
authors = [
    {name = "BiobitWorks", email = "contact@biobitworks.com"}
]
keywords = ["ai", "research", "integrity", "immune-system", "llm"]
classifiers = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Artificial Intelligence",
]
dependencies = [
    "flask>=3.0.0",
    "flask-socketio>=5.0.0",
    "requests>=2.28.0",
    "click>=8.0.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
    "mypy>=1.0.0",
]
ollama = [
    "ollama>=0.1.0",
]

[project.scripts]
immunos = "immunos.cli:main"
immunos-dashboard = "immunos.dashboard:main"
immunos-recover = "immunos.recover:main"

[project.urls]
Homepage = "https://github.com/biobitworks/immunos"
Documentation = "https://immunos.readthedocs.io"
Repository = "https://github.com/biobitworks/immunos.git"
Changelog = "https://github.com/biobitworks/immunos/blob/main/CHANGELOG.md"
Issues = "https://github.com/biobitworks/immunos/issues"

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools.package-data]
immunos = [
    "templates/*.html",
    "static/**/*",
    "kb/system/*.md",
    "config/*.json",
]
```

### MANIFEST.in

```
include LICENSE
include README.md
include CHANGELOG.md
recursive-include src/immunos/templates *
recursive-include src/immunos/static *
recursive-include src/immunos/kb *
recursive-include src/immunos/config *
```

### Version Management

```python
# src/immunos/__init__.py
__version__ = "0.1.0"
__author__ = "BiobitWorks"
```

### Publishing

```bash
# Build
python -m build

# Test upload (TestPyPI)
python -m twine upload --repository testpypi dist/*

# Production upload
python -m twine upload dist/*
```

---

## Homebrew

### Formula Structure

```ruby
# Formula/immunos.rb
class Immunos < Formula
  include Language::Python::Virtualenv

  desc "AI-assisted research integrity system"
  homepage "https://github.com/biobitworks/immunos"
  url "https://files.pythonhosted.org/packages/.../immunos-0.1.0.tar.gz"
  sha256 "abc123..."
  license "MIT"
  head "https://github.com/biobitworks/immunos.git", branch: "main"

  depends_on "python@3.11"

  resource "flask" do
    url "https://files.pythonhosted.org/packages/.../Flask-3.0.0.tar.gz"
    sha256 "..."
  end

  # Additional resources for dependencies...

  def install
    virtualenv_install_with_resources
  end

  def caveats
    <<~EOS
      IMMUNOS user data is stored in ~/.immunos/

      To start the dashboard:
        immunos-dashboard --port 5001

      To recover context:
        immunos-recover
    EOS
  end

  test do
    system "#{bin}/immunos", "--version"
  end
end
```

### Tap Repository

```
homebrew-immunos/
├── Formula/
│   └── immunos.rb
├── README.md
└── LICENSE
```

### Installation

```bash
# Add tap
brew tap biobitworks/immunos

# Install
brew install immunos

# Or direct
brew install biobitworks/immunos/immunos
```

---

## npm (Node.js)

For JavaScript components (dashboard frontend, CLI wrapper):

### package.json

```json
{
  "name": "@biobitworks/immunos",
  "version": "0.1.0",
  "description": "IMMUNOS Dashboard - AI-assisted research integrity",
  "main": "dist/index.js",
  "types": "dist/index.d.ts",
  "bin": {
    "immunos-ui": "./bin/immunos-ui.js"
  },
  "scripts": {
    "build": "tsc",
    "test": "jest",
    "lint": "eslint src/",
    "prepublishOnly": "npm run build"
  },
  "repository": {
    "type": "git",
    "url": "https://github.com/biobitworks/immunos.git"
  },
  "keywords": [
    "ai",
    "research",
    "dashboard",
    "immunos"
  ],
  "author": "BiobitWorks",
  "license": "MIT",
  "bugs": {
    "url": "https://github.com/biobitworks/immunos/issues"
  },
  "homepage": "https://github.com/biobitworks/immunos#readme",
  "engines": {
    "node": ">=18.0.0"
  },
  "files": [
    "dist/",
    "bin/",
    "static/",
    "templates/"
  ],
  "dependencies": {
    "express": "^4.18.0",
    "socket.io": "^4.6.0"
  },
  "devDependencies": {
    "@types/node": "^20.0.0",
    "typescript": "^5.0.0",
    "jest": "^29.0.0",
    "eslint": "^8.0.0"
  }
}
```

### Publishing

```bash
# Login
npm login

# Publish
npm publish --access public
```

---

## APT (Debian/Ubuntu)

### Directory Structure

```
immunos-0.1.0/
├── debian/
│   ├── control           # Package metadata
│   ├── rules             # Build rules
│   ├── changelog         # Debian changelog
│   ├── copyright         # License info
│   ├── compat            # Debhelper compat
│   ├── install           # File installation
│   └── postinst          # Post-install script
├── src/
└── ...
```

### debian/control

```
Source: immunos
Section: science
Priority: optional
Maintainer: BiobitWorks <contact@biobitworks.com>
Build-Depends: debhelper-compat (= 13),
               dh-python,
               python3-all,
               python3-setuptools
Standards-Version: 4.6.0
Homepage: https://github.com/biobitworks/immunos
Vcs-Browser: https://github.com/biobitworks/immunos
Vcs-Git: https://github.com/biobitworks/immunos.git
Rules-Requires-Root: no

Package: immunos
Architecture: all
Depends: ${python3:Depends},
         ${misc:Depends},
         python3-flask,
         python3-requests
Description: AI-assisted research integrity system
 IMMUNOS uses biological immune system metaphors to provide
 AI-assisted research integrity tools including context
 management, claim verification, and anomaly detection.
```

### debian/rules

```makefile
#!/usr/bin/make -f
export PYBUILD_NAME=immunos

%:
	dh $@ --with python3 --buildsystem=pybuild
```

### Building

```bash
# Build package
dpkg-buildpackage -us -uc

# Install
sudo dpkg -i ../immunos_0.1.0-1_all.deb
```

---

## RPM (Fedora/RHEL)

### Spec File

```spec
# immunos.spec
Name:           immunos
Version:        0.1.0
Release:        1%{?dist}
Summary:        AI-assisted research integrity system

License:        MIT
URL:            https://github.com/biobitworks/immunos
Source0:        %{pypi_source immunos}

BuildArch:      noarch
BuildRequires:  python3-devel
BuildRequires:  python3-setuptools
Requires:       python3-flask
Requires:       python3-requests

%description
IMMUNOS uses biological immune system metaphors to provide
AI-assisted research integrity tools.

%prep
%autosetup -n immunos-%{version}

%build
%py3_build

%install
%py3_install

%files
%license LICENSE
%doc README.md
%{_bindir}/immunos
%{_bindir}/immunos-dashboard
%{python3_sitelib}/immunos/
%{python3_sitelib}/immunos-%{version}*.egg-info/

%changelog
* Wed Dec 25 2024 BiobitWorks <contact@biobitworks.com> - 0.1.0-1
- Initial package
```

### Building

```bash
# Build SRPM
rpmbuild -bs immunos.spec

# Build RPM
rpmbuild -bb immunos.spec
```

---

## Docker

### Dockerfile

```dockerfile
FROM python:3.11-slim

LABEL maintainer="BiobitWorks <contact@biobitworks.com>"
LABEL version="0.1.0"
LABEL description="IMMUNOS - AI-assisted research integrity"

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application
COPY src/immunos/ ./immunos/
COPY templates/ ./templates/
COPY static/ ./static/

# Create user data directory
RUN mkdir -p /data/.immunos

ENV IMMUNOS_USER_DATA=/data/.immunos
ENV IMMUNOS_HOME=/app

EXPOSE 5001

VOLUME ["/data"]

ENTRYPOINT ["python", "-m", "immunos.dashboard"]
CMD ["--port", "5001", "--host", "0.0.0.0"]
```

### docker-compose.yml

```yaml
version: '3.8'

services:
  immunos:
    build: .
    ports:
      - "5001:5001"
    volumes:
      - immunos-data:/data
      - ./projects:/projects:ro
    environment:
      - IMMUNOS_USER_DATA=/data/.immunos
      - OLLAMA_HOST=ollama:11434
    depends_on:
      - ollama

  ollama:
    image: ollama/ollama:latest
    ports:
      - "11434:11434"
    volumes:
      - ollama-models:/root/.ollama

volumes:
  immunos-data:
  ollama-models:
```

---

## Common Requirements

### All Packages Must Include

| Item | Description |
|------|-------------|
| LICENSE | Full license text (MIT recommended) |
| README.md | Description, installation, usage |
| CHANGELOG.md | Version history (Keep a Changelog format) |
| Version string | Semantic versioning (MAJOR.MINOR.PATCH) |

### Semantic Versioning

```
MAJOR.MINOR.PATCH[-PRERELEASE][+BUILD]

Examples:
  0.1.0        - Initial development
  0.1.0-alpha  - Alpha release
  0.1.0-beta.1 - First beta
  1.0.0        - First stable release
  1.0.1        - Patch release
  1.1.0        - Minor feature release
  2.0.0        - Breaking changes
```

### CHANGELOG.md Format

```markdown
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- New feature X

### Changed
- Modified behavior Y

### Fixed
- Bug fix Z

## [0.1.0] - 2025-12-25

### Added
- Initial release
- Dashboard with T Cell Memory panel
- Context recovery system
- Snapshot management
```

---

## CI/CD Integration

### GitHub Actions (.github/workflows/publish.yml)

```yaml
name: Publish Package

on:
  release:
    types: [published]

jobs:
  pypi:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - name: Build
        run: |
          pip install build twine
          python -m build
      - name: Publish to PyPI
        env:
          TWINE_USERNAME: __token__
          TWINE_PASSWORD: ${{ secrets.PYPI_TOKEN }}
        run: twine upload dist/*

  homebrew:
    needs: pypi
    runs-on: ubuntu-latest
    steps:
      - name: Update Homebrew formula
        run: |
          # Script to update formula with new version/sha256

  docker:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: docker/build-push-action@v5
        with:
          push: true
          tags: biobitworks/immunos:${{ github.ref_name }}
```

---

## License Compatibility

| License | PyPI | npm | Homebrew | apt | rpm |
|---------|------|-----|----------|-----|-----|
| MIT | ✓ | ✓ | ✓ | ✓ | ✓ |
| Apache-2.0 | ✓ | ✓ | ✓ | ✓ | ✓ |
| GPL-3.0 | ✓ | ✓ | ✓ | ✓ | ✓ |
| BSD-3 | ✓ | ✓ | ✓ | ✓ | ✓ |

**Recommendation**: MIT License for maximum compatibility.

---

## Operating System Support

### macOS

| Version | Architectures | Package Managers |
|---------|---------------|------------------|
| macOS 12+ (Monterey) | arm64 (Apple Silicon), x86_64 (Intel) | Homebrew, pip |
| macOS 11 (Big Sur) | arm64, x86_64 | Homebrew, pip |
| macOS 10.15 (Catalina) | x86_64 only | Homebrew, pip |

#### Apple Silicon (M1/M2/M3/M4)

```bash
# Check architecture
uname -m  # Returns: arm64

# Homebrew location
/opt/homebrew/bin/brew

# Python location
/opt/homebrew/bin/python3

# Rosetta 2 for Intel binaries
softwareupdate --install-rosetta
```

#### Intel Mac

```bash
# Check architecture
uname -m  # Returns: x86_64

# Homebrew location
/usr/local/bin/brew

# Python location
/usr/local/bin/python3
```

#### Universal Binary (Fat Binary)

For compiled components, build universal binaries:
```bash
# Build for both architectures
python setup.py bdist_wheel --plat-name macosx_11_0_universal2
```

### Windows

| Version | Architectures | Package Managers |
|---------|---------------|------------------|
| Windows 11 | x64, ARM64 | pip, winget, chocolatey, scoop |
| Windows 10 | x64, x86 | pip, chocolatey, scoop |

#### winget (Windows Package Manager)

```yaml
# manifests/b/BiobitWorks/Immunos/0.1.0/BiobitWorks.Immunos.yaml
PackageIdentifier: BiobitWorks.Immunos
PackageVersion: 0.1.0
PackageLocale: en-US
Publisher: BiobitWorks
PackageName: IMMUNOS
License: MIT
ShortDescription: AI-assisted research integrity system
Installers:
  - Architecture: x64
    InstallerType: exe
    InstallerUrl: https://github.com/biobitworks/immunos/releases/download/v0.1.0/immunos-0.1.0-win64.exe
    InstallerSha256: abc123...
  - Architecture: arm64
    InstallerType: exe
    InstallerUrl: https://github.com/biobitworks/immunos/releases/download/v0.1.0/immunos-0.1.0-win-arm64.exe
    InstallerSha256: def456...
ManifestType: singleton
ManifestVersion: 1.4.0
```

#### Chocolatey

```powershell
# Install
choco install immunos

# Package spec (immunos.nuspec)
```

```xml
<?xml version="1.0" encoding="utf-8"?>
<package xmlns="http://schemas.microsoft.com/packaging/2015/06/nuspec.xsd">
  <metadata>
    <id>immunos</id>
    <version>0.1.0</version>
    <title>IMMUNOS</title>
    <authors>BiobitWorks</authors>
    <projectUrl>https://github.com/biobitworks/immunos</projectUrl>
    <licenseUrl>https://github.com/biobitworks/immunos/blob/main/LICENSE</licenseUrl>
    <requireLicenseAcceptance>false</requireLicenseAcceptance>
    <description>AI-assisted research integrity system</description>
    <tags>ai research python dashboard</tags>
  </metadata>
</package>
```

#### Windows Paths

```
%LOCALAPPDATA%\immunos\           # User data (C:\Users\X\AppData\Local\immunos)
%PROGRAMFILES%\immunos\           # System installation (C:\Program Files\immunos)
```

### Linux

| Distribution | Package Format | Package Managers |
|--------------|----------------|------------------|
| Ubuntu/Debian | .deb | apt, pip, snap |
| Fedora/RHEL | .rpm | dnf/yum, pip |
| Arch Linux | PKGBUILD | pacman, pip, AUR |
| Alpine | .apk | apk, pip |

#### Snap Package

```yaml
# snapcraft.yaml
name: immunos
base: core22
version: '0.1.0'
summary: AI-assisted research integrity system
description: |
  IMMUNOS uses biological immune system metaphors for
  AI-assisted research integrity.
grade: stable
confinement: strict

apps:
  immunos:
    command: bin/immunos
    plugs:
      - network
      - home
  dashboard:
    command: bin/immunos-dashboard
    plugs:
      - network
      - network-bind

parts:
  immunos:
    plugin: python
    source: .
    python-packages:
      - flask
      - flask-socketio
```

#### Flatpak

```yaml
# com.biobitworks.Immunos.yaml
app-id: com.biobitworks.Immunos
runtime: org.freedesktop.Platform
runtime-version: '23.08'
sdk: org.freedesktop.Sdk
command: immunos-dashboard

finish-args:
  - --share=network
  - --socket=wayland
  - --socket=fallback-x11

modules:
  - name: immunos
    buildsystem: simple
    build-commands:
      - pip3 install --prefix=/app .
    sources:
      - type: git
        url: https://github.com/biobitworks/immunos.git
        tag: v0.1.0
```

#### AppImage

```yaml
# AppImageBuilder.yml
version: 1
AppDir:
  path: ./AppDir
  app_info:
    id: com.biobitworks.Immunos
    name: IMMUNOS
    icon: immunos
    version: 0.1.0
    exec: usr/bin/python3
    exec_args: "$APPDIR/usr/bin/immunos-dashboard $@"

  apt:
    arch: amd64
    sources:
      - sourceline: deb http://archive.ubuntu.com/ubuntu jammy main
    include:
      - python3
      - python3-pip
```

---

## Processor Architecture Support

### Architecture Matrix

| Architecture | Identifier | Platforms | Notes |
|--------------|------------|-----------|-------|
| **x86_64** | amd64, x64 | Intel/AMD 64-bit | Most common desktop/server |
| **ARM64** | aarch64, arm64 | Apple Silicon, AWS Graviton, Raspberry Pi 4+ | Growing server adoption |
| **x86** | i386, i686 | Intel/AMD 32-bit | Legacy, declining support |
| **ARMv7** | armhf | Raspberry Pi 3, embedded | IoT/embedded only |

### Python Wheel Tags

```
immunos-0.1.0-py3-none-any.whl           # Pure Python (all platforms)
immunos-0.1.0-cp311-cp311-manylinux_2_17_x86_64.whl    # Linux x86_64
immunos-0.1.0-cp311-cp311-manylinux_2_17_aarch64.whl   # Linux ARM64
immunos-0.1.0-cp311-cp311-macosx_11_0_arm64.whl        # macOS ARM64
immunos-0.1.0-cp311-cp311-macosx_10_9_x86_64.whl       # macOS Intel
immunos-0.1.0-cp311-cp311-macosx_11_0_universal2.whl   # macOS Universal
immunos-0.1.0-cp311-cp311-win_amd64.whl                # Windows x64
immunos-0.1.0-cp311-cp311-win_arm64.whl                # Windows ARM64
```

### Docker Multi-Architecture

```dockerfile
# Build for multiple architectures
FROM --platform=$TARGETPLATFORM python:3.11-slim

ARG TARGETPLATFORM
ARG BUILDPLATFORM

RUN echo "Building on $BUILDPLATFORM for $TARGETPLATFORM"
```

```bash
# Build and push multi-arch image
docker buildx create --use
docker buildx build --platform linux/amd64,linux/arm64 \
  -t biobitworks/immunos:0.1.0 --push .
```

### Ollama Architecture Support

| Architecture | Ollama Support | GPU Support |
|--------------|----------------|-------------|
| x86_64 Linux | ✓ Full | CUDA, ROCm |
| ARM64 Linux | ✓ Full | Limited |
| x86_64 macOS | ✓ Full | Metal (limited) |
| ARM64 macOS | ✓ Full | Metal (native) |
| x86_64 Windows | ✓ Full | CUDA |
| ARM64 Windows | ✓ Preview | DirectML |

### Build Matrix (GitHub Actions)

```yaml
# .github/workflows/build.yml
name: Build

on: [push, pull_request]

jobs:
  build:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python-version: ['3.9', '3.10', '3.11', '3.12']
        include:
          - os: macos-latest
            arch: arm64
          - os: macos-13  # Intel runner
            arch: x86_64
          - os: ubuntu-latest
            arch: x86_64
          # ARM64 Linux via QEMU
          - os: ubuntu-latest
            arch: arm64
            platform: linux/arm64

    runs-on: ${{ matrix.os }}

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Set up QEMU (for ARM64)
        if: matrix.arch == 'arm64' && runner.os == 'Linux'
        uses: docker/setup-qemu-action@v3

      - name: Build wheel
        run: |
          pip install build
          python -m build

      - name: Upload artifact
        uses: actions/upload-artifact@v4
        with:
          name: wheel-${{ matrix.os }}-${{ matrix.python-version }}-${{ matrix.arch }}
          path: dist/*.whl
```

---

## Platform-Specific Paths

### Default Installation Paths

| Platform | System Install | User Data |
|----------|---------------|-----------|
| **macOS (Homebrew)** | `/opt/homebrew/lib/python3.x/site-packages/` | `~/.immunos/` |
| **macOS (pip user)** | `~/Library/Python/3.x/lib/python/site-packages/` | `~/.immunos/` |
| **Linux (pip)** | `/usr/local/lib/python3.x/dist-packages/` | `~/.immunos/` |
| **Linux (apt)** | `/usr/lib/python3/dist-packages/` | `~/.immunos/` |
| **Windows (pip)** | `%PROGRAMFILES%\Python3x\Lib\site-packages\` | `%LOCALAPPDATA%\immunos\` |
| **Windows (user)** | `%APPDATA%\Python\Python3x\site-packages\` | `%LOCALAPPDATA%\immunos\` |

### XDG Base Directory (Linux)

```python
import os
from pathlib import Path

def get_user_data_dir():
    """Get user data directory following XDG spec."""
    if os.name == 'nt':  # Windows
        return Path(os.environ.get('LOCALAPPDATA', Path.home() / 'AppData' / 'Local')) / 'immunos'
    elif sys.platform == 'darwin':  # macOS
        return Path.home() / '.immunos'
    else:  # Linux/Unix
        xdg_data = os.environ.get('XDG_DATA_HOME', Path.home() / '.local' / 'share')
        return Path(xdg_data) / 'immunos'

def get_user_config_dir():
    """Get user config directory following XDG spec."""
    if os.name == 'nt':
        return Path(os.environ.get('APPDATA', Path.home() / 'AppData' / 'Roaming')) / 'immunos'
    elif sys.platform == 'darwin':
        return Path.home() / '.immunos' / 'config'
    else:
        xdg_config = os.environ.get('XDG_CONFIG_HOME', Path.home() / '.config')
        return Path(xdg_config) / 'immunos'

def get_user_cache_dir():
    """Get user cache directory."""
    if os.name == 'nt':
        return Path(os.environ.get('LOCALAPPDATA', '')) / 'immunos' / 'cache'
    elif sys.platform == 'darwin':
        return Path.home() / 'Library' / 'Caches' / 'immunos'
    else:
        xdg_cache = os.environ.get('XDG_CACHE_HOME', Path.home() / '.cache')
        return Path(xdg_cache) / 'immunos'
```

---

## Testing Matrix

| OS | Arch | Python | Priority |
|----|------|--------|----------|
| Ubuntu 22.04 | x86_64 | 3.11 | P0 (primary) |
| macOS 14 | arm64 | 3.11 | P0 |
| Windows 11 | x64 | 3.11 | P0 |
| Ubuntu 22.04 | arm64 | 3.11 | P1 |
| macOS 13 | x86_64 | 3.11 | P1 |
| Fedora 39 | x86_64 | 3.12 | P2 |
| Debian 12 | x86_64 | 3.11 | P2 |
| Windows 11 | arm64 | 3.11 | P2 |
| Alpine 3.19 | x86_64 | 3.11 | P3 (Docker) |

---

**Last Updated**: 2025-12-25
**Status**: Planning
**Target**: v0.1.0 alpha release
