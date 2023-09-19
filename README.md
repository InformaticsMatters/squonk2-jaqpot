# squonk2-jaqpot

![Data Manager Job: 2021.1](https://img.shields.io/badge/data%20manager%20job-2021.1-000000?labelColor=dc332e)

[![build](https://github.com/InformaticsMatters/squonk2-jaqpot/actions/workflows/build.yaml/badge.svg)](https://github.com/InformaticsMatters/squonk2-jaqpot/actions/workflows/build.yaml)
[![publish-tag](https://github.com/InformaticsMatters/squonk2-jaqpot/actions/workflows/publish-tag.yaml/badge.svg)](https://github.com/InformaticsMatters/squonk2-jaqpot/actions/workflows/publish-tag.yaml)
[![publish-stable](https://github.com/InformaticsMatters/squonk2-jaqpot/actions/workflows/publish-stable.yaml/badge.svg)](https://github.com/InformaticsMatters/squonk2-jaqpot/actions/workflows/publish-stable.yaml)

![GitHub](https://img.shields.io/github/license/informaticsmatters/squonk2-jaqpot)

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/informaticsmatters/squonk2-jaqpot)

A Squonk2 Data Manager Job repository for Jaqpot (Torch) models.

From a fork you should be able to build and run the tests for the example
Job that it defines, start with this, and you'll know you're starting with
a working framework: -

    python -m venv venv
    source venv/bin/activate
    python -m pip install -r build-requirements.txt

    docker-compose build
    jote

    deactivate

---

[buildx]: https://docs.docker.com/buildx/working-with-buildx
[buildx gist]: https://gist.github.com/alanbchristie/14da3444f3fed6f0adcf877a82b56804.js
[im-jote]: https://pypi.org/project/im-jote
