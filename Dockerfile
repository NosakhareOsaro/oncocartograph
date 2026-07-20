# syntax=docker/dockerfile:1
#
# NOTE: this image currently covers only the Python side of the pipeline
# (package install, lint/type-check/test). The MOFA+ integration stage
# requires an R environment (R + MOFA2) and will extend this image (or add
# a sibling image/service in docker-compose.yml) once feat/mofa-integration
# lands -- it is intentionally not stubbed out here in advance.
FROM python:3.11-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src ./src
COPY tests ./tests

RUN pip install --upgrade pip \
    && pip install -e ".[dev]"

CMD ["pytest"]
